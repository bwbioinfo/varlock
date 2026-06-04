use std::{
    collections::{HashMap, HashSet},
    fs::File,
    io::{BufRead, Write},
    path::Path,
    time::Instant,
};

use anyhow::{Context, Result, bail};
use noodles_bgzf as bgzf;
use noodles_csi::binning_index::index::reference_sequence::bin::Chunk;

use crate::{
    ExecutionContext, FilterArgs, IndexType, log_verbose,
    vcf::{
        IndexRecord, OutputIndex, VcfRecord, open_text_reader, parse_header_samples, parse_info_ref,
    },
};

#[derive(Clone, Debug)]
struct MaxInfoFilter {
    field: String,
    max: f64,
}

#[derive(Clone, Debug)]
struct SampleGtFilter {
    sample: String,
    gt: String,
}

#[derive(Clone, Debug)]
struct SampleMinDpFilter {
    sample: String,
    min_dp: u32,
}

#[derive(Clone, Debug)]
struct GroupGtFilter {
    group: String,
    gt: String,
}

#[derive(Clone, Debug)]
struct GroupMinDpFilter {
    group: String,
    min_dp: u32,
}

#[derive(Debug)]
struct FilterSpec {
    require_info: HashSet<String>,
    exclude_info: HashSet<String>,
    max_info: Vec<MaxInfoFilter>,
    expressions: Vec<Expr>,
    sample_groups: HashMap<String, Vec<String>>,
    sample_has_alt: Vec<String>,
    sample_gt: Vec<SampleGtFilter>,
    sample_min_dp: Vec<SampleMinDpFilter>,
    group_any_has_alt: Vec<String>,
    group_any_gt: Vec<GroupGtFilter>,
    group_all_min_dp: Vec<GroupMinDpFilter>,
}

pub(crate) fn run(args: FilterArgs, ctx: &ExecutionContext) -> Result<()> {
    let started = Instant::now();
    let spec = FilterSpec::from_args(&args)?;
    let metrics = filter_vcf(&args.input, &args.output, args.index_type, &spec)?;

    log_verbose(
        ctx,
        format!(
            "filter stage=done input_records={} output_records={} output={} elapsed={:.2?}",
            metrics.input_records,
            metrics.output_records,
            args.output.display(),
            started.elapsed()
        ),
    );
    Ok(())
}

impl FilterSpec {
    fn from_args(args: &FilterArgs) -> Result<Self> {
        let max_info = args
            .max_info
            .iter()
            .map(|value| {
                let (field, max) = value.split_once('=').with_context(|| {
                    format!("invalid --max-info {value:?}; expected FIELD=VALUE")
                })?;
                if field.is_empty() || max.is_empty() {
                    bail!("invalid --max-info {value:?}; expected FIELD=VALUE");
                }
                let max = max
                    .parse::<f64>()
                    .with_context(|| format!("invalid --max-info threshold {max:?}"))?;
                Ok(MaxInfoFilter {
                    field: field.to_string(),
                    max,
                })
            })
            .collect::<Result<Vec<_>>>()?;
        let sample_groups = args
            .sample_groups
            .iter()
            .map(|value| {
                let (name, samples) = split_name_value(value, "--sample-group")?;
                let samples = samples
                    .split(',')
                    .filter(|sample| !sample.is_empty())
                    .map(str::to_string)
                    .collect::<Vec<_>>();
                if samples.is_empty() {
                    bail!("invalid --sample-group {value:?}; expected NAME=SAMPLE[,SAMPLE...]");
                }
                Ok((name.to_string(), samples))
            })
            .collect::<Result<HashMap<_, _>>>()?;
        let sample_gt = args
            .sample_gt
            .iter()
            .map(|value| {
                let (sample, gt) = split_name_value(value, "--sample-gt")?;
                Ok(SampleGtFilter {
                    sample: sample.to_string(),
                    gt: gt.to_string(),
                })
            })
            .collect::<Result<Vec<_>>>()?;
        let sample_min_dp = args
            .sample_min_dp
            .iter()
            .map(|value| {
                let (sample, min_dp) = split_name_value(value, "--sample-min-dp")?;
                Ok(SampleMinDpFilter {
                    sample: sample.to_string(),
                    min_dp: parse_dp(min_dp, "--sample-min-dp")?,
                })
            })
            .collect::<Result<Vec<_>>>()?;
        let group_any_gt = args
            .group_any_gt
            .iter()
            .map(|value| {
                let (group, gt) = split_name_value(value, "--group-any-gt")?;
                Ok(GroupGtFilter {
                    group: group.to_string(),
                    gt: gt.to_string(),
                })
            })
            .collect::<Result<Vec<_>>>()?;
        let group_all_min_dp = args
            .group_all_min_dp
            .iter()
            .map(|value| {
                let (group, min_dp) = split_name_value(value, "--group-all-min-dp")?;
                Ok(GroupMinDpFilter {
                    group: group.to_string(),
                    min_dp: parse_dp(min_dp, "--group-all-min-dp")?,
                })
            })
            .collect::<Result<Vec<_>>>()?;
        let expressions = args
            .expressions
            .iter()
            .map(|value| parse_expr(value))
            .collect::<Result<Vec<_>>>()?;

        Ok(Self {
            require_info: args.require_info.iter().cloned().collect(),
            exclude_info: args.exclude_info.iter().cloned().collect(),
            max_info,
            expressions,
            sample_groups,
            sample_has_alt: args.sample_has_alt.clone(),
            sample_gt,
            sample_min_dp,
            group_any_has_alt: args.group_any_has_alt.clone(),
            group_any_gt,
            group_all_min_dp,
        })
    }

    fn validate_samples(&self, sample_names: &[String]) -> Result<()> {
        let known = sample_names
            .iter()
            .map(String::as_str)
            .collect::<HashSet<_>>();
        for sample in &self.sample_has_alt {
            validate_known_sample(sample, &known)?;
        }
        for filter in &self.sample_gt {
            validate_known_sample(&filter.sample, &known)?;
        }
        for filter in &self.sample_min_dp {
            validate_known_sample(&filter.sample, &known)?;
        }
        for (group, samples) in &self.sample_groups {
            if samples.is_empty() {
                bail!("sample group {group:?} has no samples");
            }
            for sample in samples {
                validate_known_sample(sample, &known)
                    .with_context(|| format!("invalid sample group {group:?}"))?;
            }
        }
        for group in &self.group_any_has_alt {
            self.group_samples(group)?;
        }
        for filter in &self.group_any_gt {
            self.group_samples(&filter.group)?;
        }
        for filter in &self.group_all_min_dp {
            self.group_samples(&filter.group)?;
        }
        Ok(())
    }

    fn keep_record(&self, record: &VcfRecord<'_>) -> Result<bool> {
        let info = parse_info_ref(record.info_text);
        for field in &self.require_info {
            if !info.contains_key(field.as_str()) {
                return Ok(false);
            }
        }
        for field in &self.exclude_info {
            if info.contains_key(field.as_str()) {
                return Ok(false);
            }
        }
        for filter in &self.max_info {
            let Some(value) = info.get(filter.field.as_str()) else {
                return Ok(false);
            };
            if !all_numeric_values_at_most(value, filter.max)? {
                return Ok(false);
            }
        }
        for expr in &self.expressions {
            if !eval_expr(expr, &info)? {
                return Ok(false);
            }
        }
        for sample in &self.sample_has_alt {
            if !sample_has_alt(record, sample)? {
                return Ok(false);
            }
        }
        for filter in &self.sample_gt {
            if record.sample_field(&filter.sample, "GT") != Some(filter.gt.as_str()) {
                return Ok(false);
            }
        }
        for filter in &self.sample_min_dp {
            if !sample_min_dp(record, &filter.sample, filter.min_dp)? {
                return Ok(false);
            }
        }
        for group in &self.group_any_has_alt {
            let mut any = false;
            for sample in self.group_samples(group)? {
                if sample_has_alt(record, sample)? {
                    any = true;
                    break;
                }
            }
            if !any {
                return Ok(false);
            }
        }
        for filter in &self.group_any_gt {
            if !self
                .group_samples(&filter.group)?
                .iter()
                .any(|sample| record.sample_field(sample, "GT") == Some(filter.gt.as_str()))
            {
                return Ok(false);
            }
        }
        for filter in &self.group_all_min_dp {
            for sample in self.group_samples(&filter.group)? {
                if !sample_min_dp(record, sample, filter.min_dp)? {
                    return Ok(false);
                }
            }
        }
        Ok(true)
    }

    fn group_samples(&self, group: &str) -> Result<&[String]> {
        self.sample_groups
            .get(group)
            .map(Vec::as_slice)
            .with_context(|| format!("sample group {group:?} is not defined"))
    }
}

#[derive(Default)]
struct FilterMetrics {
    input_records: usize,
    output_records: usize,
}

fn filter_vcf(
    input: &Path,
    output: &Path,
    index_type: IndexType,
    spec: &FilterSpec,
) -> Result<FilterMetrics> {
    if let Some(parent) = output.parent()
        && !parent.as_os_str().is_empty()
    {
        std::fs::create_dir_all(parent)
            .with_context(|| format!("failed to create output directory {}", parent.display()))?;
    }

    let mut reader = open_text_reader(input)
        .with_context(|| format!("failed to open input VCF {}", input.display()))?;
    let output_file = File::create(output)
        .with_context(|| format!("failed to create output {}", output.display()))?;
    let mut writer = bgzf::io::writer::Builder::default().build_from_writer(output_file);
    let mut output_index = OutputIndex::new(index_type);

    let mut line = String::new();
    let mut saw_column_header = false;
    let mut sample_index = HashMap::new();
    let mut metrics = FilterMetrics::default();
    while reader.read_line(&mut line)? != 0 {
        let trimmed = line.trim_end_matches(['\r', '\n']);
        if trimmed.starts_with("#CHROM") {
            saw_column_header = true;
            let (sample_names, parsed_sample_index) = parse_header_samples(trimmed);
            sample_index = parsed_sample_index;
            spec.validate_samples(&sample_names)?;
            writeln!(writer, "{trimmed}")?;
        } else if trimmed.starts_with('#') {
            writeln!(writer, "{trimmed}")?;
        } else {
            metrics.input_records += 1;
            let fields = trimmed.split('\t').collect::<Vec<_>>();
            if fields.len() < 8 {
                bail!("invalid VCF record with fewer than 8 fields: {trimmed}");
            }
            let record = VcfRecord::new(&fields, &sample_index);
            if spec.keep_record(&record)? {
                let record = IndexRecord::from_fields(&fields)?;
                let chunk_start = writer.virtual_position();
                writeln!(writer, "{trimmed}")?;
                let chunk_end = writer.virtual_position();
                output_index.add_record(&record, Chunk::new(chunk_start, chunk_end))?;
                metrics.output_records += 1;
            }
        }
        line.clear();
    }
    if !saw_column_header {
        bail!("input VCF is missing #CHROM header line");
    }
    writer
        .try_finish()
        .context("failed to finish bgzip output")?;
    output_index.write(output)?;
    Ok(metrics)
}

fn split_name_value<'a>(value: &'a str, flag: &str) -> Result<(&'a str, &'a str)> {
    let (name, parsed_value) = value
        .split_once('=')
        .with_context(|| format!("invalid {flag} {value:?}; expected NAME=VALUE"))?;
    if name.is_empty() || parsed_value.is_empty() {
        bail!("invalid {flag} {value:?}; expected NAME=VALUE");
    }
    Ok((name, parsed_value))
}

fn parse_dp(value: &str, flag: &str) -> Result<u32> {
    value
        .parse::<u32>()
        .with_context(|| format!("invalid {flag} depth threshold {value:?}"))
}

fn validate_known_sample(sample: &str, known: &HashSet<&str>) -> Result<()> {
    if known.contains(sample) {
        Ok(())
    } else {
        bail!("sample {sample:?} was not found in VCF header")
    }
}

fn sample_has_alt(record: &VcfRecord<'_>, sample: &str) -> Result<bool> {
    let Some(gt) = record.sample_field(sample, "GT") else {
        return Ok(false);
    };
    genotype_has_alt(gt)
}

fn sample_min_dp(record: &VcfRecord<'_>, sample: &str, min_dp: u32) -> Result<bool> {
    let Some(dp) = record.sample_field(sample, "DP") else {
        return Ok(false);
    };
    let dp = dp
        .parse::<u32>()
        .with_context(|| format!("FORMAT/DP for sample {sample:?} is not an integer"))?;
    Ok(dp >= min_dp)
}

fn genotype_has_alt(gt: &str) -> Result<bool> {
    if gt == "." || gt == "./." || gt == ".|." {
        return Ok(false);
    }
    for allele in gt.split(['/', '|']) {
        if allele == "." || allele.is_empty() {
            continue;
        }
        let allele = allele
            .parse::<u32>()
            .with_context(|| format!("FORMAT/GT allele {allele:?} is not an integer"))?;
        if allele > 0 {
            return Ok(true);
        }
    }
    Ok(false)
}

#[derive(Clone, Debug)]
enum Expr {
    Field(String),
    Number(f64),
    String(String),
    Missing(String),
    Not(Box<Expr>),
    Compare {
        left: Box<Expr>,
        op: CompareOp,
        right: Box<Expr>,
    },
    And(Box<Expr>, Box<Expr>),
    Or(Box<Expr>, Box<Expr>),
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum CompareOp {
    Eq,
    Ne,
    Lt,
    Le,
    Gt,
    Ge,
}

#[derive(Clone, Debug, PartialEq)]
enum Token {
    Ident(String),
    Number(f64),
    String(String),
    Missing,
    And,
    Or,
    Not,
    Eq,
    Ne,
    Lt,
    Le,
    Gt,
    Ge,
    LParen,
    RParen,
    Comma,
}

fn parse_expr(input: &str) -> Result<Expr> {
    let tokens = tokenize_expr(input)?;
    if tokens.is_empty() {
        bail!("invalid --expr {input:?}; expression cannot be empty");
    }
    let mut parser = ExprParser { tokens, pos: 0 };
    let expr = parser.parse_or()?;
    if parser.peek().is_some() {
        bail!("invalid --expr {input:?}; unexpected token after expression");
    }
    Ok(expr)
}

struct ExprParser {
    tokens: Vec<Token>,
    pos: usize,
}

impl ExprParser {
    fn parse_or(&mut self) -> Result<Expr> {
        let mut expr = self.parse_and()?;
        while self.consume(|token| matches!(token, Token::Or)).is_some() {
            let right = self.parse_and()?;
            expr = Expr::Or(Box::new(expr), Box::new(right));
        }
        Ok(expr)
    }

    fn parse_and(&mut self) -> Result<Expr> {
        let mut expr = self.parse_compare()?;
        while self.consume(|token| matches!(token, Token::And)).is_some() {
            let right = self.parse_compare()?;
            expr = Expr::And(Box::new(expr), Box::new(right));
        }
        Ok(expr)
    }

    fn parse_compare(&mut self) -> Result<Expr> {
        let left = self.parse_unary()?;
        let Some(op) = self.consume_compare_op() else {
            return Ok(left);
        };
        let right = self.parse_unary()?;
        Ok(Expr::Compare {
            left: Box::new(left),
            op,
            right: Box::new(right),
        })
    }

    fn parse_unary(&mut self) -> Result<Expr> {
        if self.consume(|token| matches!(token, Token::Not)).is_some() {
            return Ok(Expr::Not(Box::new(self.parse_unary()?)));
        }
        self.parse_primary()
    }

    fn parse_primary(&mut self) -> Result<Expr> {
        match self.next() {
            Some(Token::Ident(value)) => Ok(Expr::Field(value)),
            Some(Token::Number(value)) => Ok(Expr::Number(value)),
            Some(Token::String(value)) => Ok(Expr::String(value)),
            Some(Token::Missing) => {
                self.expect_lparen("missing")?;
                let field = match self.next() {
                    Some(Token::Ident(field)) => field,
                    other => bail!("missing() expects a field name, got {other:?}"),
                };
                self.expect_rparen("missing")?;
                Ok(Expr::Missing(field))
            }
            Some(Token::LParen) => {
                let expr = self.parse_or()?;
                self.expect_rparen("parenthesized expression")?;
                Ok(expr)
            }
            other => bail!("unexpected expression token {other:?}"),
        }
    }

    fn consume_compare_op(&mut self) -> Option<CompareOp> {
        let op = match self.peek()? {
            Token::Eq => CompareOp::Eq,
            Token::Ne => CompareOp::Ne,
            Token::Lt => CompareOp::Lt,
            Token::Le => CompareOp::Le,
            Token::Gt => CompareOp::Gt,
            Token::Ge => CompareOp::Ge,
            _ => return None,
        };
        self.pos += 1;
        Some(op)
    }

    fn expect_lparen(&mut self, context: &str) -> Result<()> {
        if self
            .consume(|token| matches!(token, Token::LParen))
            .is_some()
        {
            Ok(())
        } else {
            bail!("{context} expects '('")
        }
    }

    fn expect_rparen(&mut self, context: &str) -> Result<()> {
        if self
            .consume(|token| matches!(token, Token::RParen))
            .is_some()
        {
            Ok(())
        } else {
            bail!("{context} expects ')'")
        }
    }

    fn consume(&mut self, predicate: impl FnOnce(&Token) -> bool) -> Option<Token> {
        let token = self.peek()?.clone();
        if predicate(&token) {
            self.pos += 1;
            Some(token)
        } else {
            None
        }
    }

    fn next(&mut self) -> Option<Token> {
        let token = self.peek()?.clone();
        self.pos += 1;
        Some(token)
    }

    fn peek(&self) -> Option<&Token> {
        self.tokens.get(self.pos)
    }
}

fn tokenize_expr(input: &str) -> Result<Vec<Token>> {
    let mut tokens = Vec::new();
    let mut chars = input.char_indices().peekable();
    while let Some((start, ch)) = chars.next() {
        match ch {
            c if c.is_whitespace() => {}
            '(' => tokens.push(Token::LParen),
            ')' => tokens.push(Token::RParen),
            ',' => tokens.push(Token::Comma),
            '!' => {
                if consume_char(&mut chars, '=') {
                    tokens.push(Token::Ne);
                } else {
                    tokens.push(Token::Not);
                }
            }
            '=' => {
                if consume_char(&mut chars, '=') {
                    tokens.push(Token::Eq);
                } else {
                    bail!("invalid expression operator at byte {start}; use ==");
                }
            }
            '<' => {
                if consume_char(&mut chars, '=') {
                    tokens.push(Token::Le);
                } else {
                    tokens.push(Token::Lt);
                }
            }
            '>' => {
                if consume_char(&mut chars, '=') {
                    tokens.push(Token::Ge);
                } else {
                    tokens.push(Token::Gt);
                }
            }
            '&' => {
                if consume_char(&mut chars, '&') {
                    tokens.push(Token::And);
                } else {
                    bail!("invalid expression operator at byte {start}; use &&");
                }
            }
            '|' => {
                if consume_char(&mut chars, '|') {
                    tokens.push(Token::Or);
                } else {
                    bail!("invalid expression operator at byte {start}; use ||");
                }
            }
            '"' | '\'' => tokens.push(Token::String(read_quoted(input, &mut chars, ch)?)),
            c if c.is_ascii_digit() || c == '.' => {
                let mut end = start + ch.len_utf8();
                while let Some(&(idx, next)) = chars.peek() {
                    if next.is_ascii_digit()
                        || matches!(next, '.' | 'e' | 'E' | '+' | '-')
                        || next == '_'
                    {
                        chars.next();
                        end = idx + next.len_utf8();
                    } else {
                        break;
                    }
                }
                let raw = input[start..end].replace('_', "");
                let value = raw
                    .parse::<f64>()
                    .with_context(|| format!("invalid expression number {raw:?}"))?;
                tokens.push(Token::Number(value));
            }
            c if is_ident_start(c) => {
                let mut end = start + ch.len_utf8();
                while let Some(&(idx, next)) = chars.peek() {
                    if is_ident_continue(next) {
                        chars.next();
                        end = idx + next.len_utf8();
                    } else {
                        break;
                    }
                }
                let ident = &input[start..end];
                if ident == "missing" {
                    tokens.push(Token::Missing);
                } else {
                    tokens.push(Token::Ident(ident.to_string()));
                }
            }
            _ => bail!("invalid expression character {ch:?} at byte {start}"),
        }
    }
    Ok(tokens)
}

fn consume_char(
    chars: &mut std::iter::Peekable<std::str::CharIndices<'_>>,
    expected: char,
) -> bool {
    if chars.peek().map(|(_, ch)| *ch == expected).unwrap_or(false) {
        chars.next();
        true
    } else {
        false
    }
}

fn read_quoted(
    input: &str,
    chars: &mut std::iter::Peekable<std::str::CharIndices<'_>>,
    quote: char,
) -> Result<String> {
    let mut out = String::new();
    while let Some((_, ch)) = chars.next() {
        if ch == quote {
            return Ok(out);
        }
        if ch == '\\' {
            let Some((_, escaped)) = chars.next() else {
                bail!("unterminated escape in expression string {input:?}");
            };
            out.push(escaped);
        } else {
            out.push(ch);
        }
    }
    bail!("unterminated expression string {input:?}")
}

fn is_ident_start(ch: char) -> bool {
    ch.is_ascii_alphabetic() || ch == '_'
}

fn is_ident_continue(ch: char) -> bool {
    ch.is_ascii_alphanumeric() || matches!(ch, '_' | '.' | '-')
}

fn eval_expr(expr: &Expr, info: &HashMap<&str, &str>) -> Result<bool> {
    match expr {
        Expr::Field(field) => Ok(info.contains_key(field.as_str())),
        Expr::Missing(field) => Ok(!info.contains_key(field.as_str())),
        Expr::Number(value) => Ok(*value != 0.0),
        Expr::String(value) => Ok(!value.is_empty()),
        Expr::Not(inner) => Ok(!eval_expr(inner, info)?),
        Expr::And(left, right) => Ok(eval_expr(left, info)? && eval_expr(right, info)?),
        Expr::Or(left, right) => Ok(eval_expr(left, info)? || eval_expr(right, info)?),
        Expr::Compare { left, op, right } => eval_compare(left, *op, right, info),
    }
}

fn eval_compare(
    left: &Expr,
    op: CompareOp,
    right: &Expr,
    info: &HashMap<&str, &str>,
) -> Result<bool> {
    let left = expr_value(left, info)?;
    let right = expr_value(right, info)?;
    match (left.as_ref(), right.as_ref()) {
        (Value::Missing, _) | (_, Value::Missing) => Ok(matches!(op, CompareOp::Ne)),
        (Value::Number(left), Value::Number(right)) => compare_numbers(*left, op, *right),
        (Value::String(left), Value::String(right)) => compare_strings(left, op, right),
        (Value::Number(left), Value::String(right)) => {
            if let Ok(right) = right.parse::<f64>() {
                compare_numbers(*left, op, right)
            } else {
                compare_strings(&left.to_string(), op, right)
            }
        }
        (Value::String(left), Value::Number(right)) => {
            if let Ok(left) = left.parse::<f64>() {
                compare_numbers(left, op, *right)
            } else {
                compare_strings(left, op, &right.to_string())
            }
        }
    }
}

#[derive(Clone, Debug)]
enum Value {
    Missing,
    Number(f64),
    String(String),
}

impl Value {
    fn as_ref(&self) -> &Self {
        self
    }
}

fn expr_value(expr: &Expr, info: &HashMap<&str, &str>) -> Result<Value> {
    match expr {
        Expr::Field(field) => {
            let Some(value) = info.get(field.as_str()) else {
                return Ok(Value::Missing);
            };
            if let Ok(number) = value.parse::<f64>() {
                Ok(Value::Number(number))
            } else {
                Ok(Value::String((*value).to_string()))
            }
        }
        Expr::Number(value) => Ok(Value::Number(*value)),
        Expr::String(value) => Ok(Value::String(value.clone())),
        Expr::Missing(field) => Ok(Value::Number(if info.contains_key(field.as_str()) {
            0.0
        } else {
            1.0
        })),
        Expr::Not(_) | Expr::Compare { .. } | Expr::And(_, _) | Expr::Or(_, _) => {
            Ok(Value::Number(if eval_expr(expr, info)? {
                1.0
            } else {
                0.0
            }))
        }
    }
}

fn compare_numbers(left: f64, op: CompareOp, right: f64) -> Result<bool> {
    Ok(match op {
        CompareOp::Eq => left == right,
        CompareOp::Ne => left != right,
        CompareOp::Lt => left < right,
        CompareOp::Le => left <= right,
        CompareOp::Gt => left > right,
        CompareOp::Ge => left >= right,
    })
}

fn compare_strings(left: &str, op: CompareOp, right: &str) -> Result<bool> {
    Ok(match op {
        CompareOp::Eq => left == right,
        CompareOp::Ne => left != right,
        CompareOp::Lt => left < right,
        CompareOp::Le => left <= right,
        CompareOp::Gt => left > right,
        CompareOp::Ge => left >= right,
    })
}

fn all_numeric_values_at_most(value: &str, max: f64) -> Result<bool> {
    let mut saw_value = false;
    for part in value.split(',') {
        if part == "." || part.is_empty() {
            continue;
        }
        saw_value = true;
        let parsed = part
            .parse::<f64>()
            .with_context(|| format!("INFO value {part:?} is not numeric"))?;
        if parsed > max {
            return Ok(false);
        }
    }
    Ok(saw_value)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Read;
    use tempfile::tempdir;

    fn empty_spec() -> FilterSpec {
        FilterSpec {
            require_info: HashSet::new(),
            exclude_info: HashSet::new(),
            max_info: Vec::new(),
            expressions: Vec::new(),
            sample_groups: HashMap::new(),
            sample_has_alt: Vec::new(),
            sample_gt: Vec::new(),
            sample_min_dp: Vec::new(),
            group_any_has_alt: Vec::new(),
            group_any_gt: Vec::new(),
            group_all_min_dp: Vec::new(),
        }
    }

    fn record<'a>(
        fields: &'a [&'a str],
        sample_index: &'a HashMap<String, usize>,
    ) -> VcfRecord<'a> {
        VcfRecord::new(fields, sample_index)
    }

    fn sample_index(samples: &[&str]) -> HashMap<String, usize> {
        samples
            .iter()
            .enumerate()
            .map(|(i, sample)| ((*sample).to_string(), i))
            .collect()
    }

    #[test]
    fn filter_spec_requires_and_excludes_info() -> Result<()> {
        let spec = FilterSpec {
            require_info: ["AF".to_string()].into_iter().collect(),
            exclude_info: ["COMMON".to_string()].into_iter().collect(),
            ..empty_spec()
        };
        let sample_index = HashMap::new();
        assert!(spec.keep_record(&record(
            &["chr1", "1", ".", "A", "C", ".", "PASS", "AF=0.1"],
            &sample_index
        ))?);
        assert!(!spec.keep_record(&record(
            &["chr1", "1", ".", "A", "C", ".", "PASS", "DP=10"],
            &sample_index
        ))?);
        assert!(!spec.keep_record(&record(
            &["chr1", "1", ".", "A", "C", ".", "PASS", "AF=0.1;COMMON"],
            &sample_index
        ))?);
        Ok(())
    }

    #[test]
    fn filter_spec_applies_max_info_to_all_numeric_values() -> Result<()> {
        let spec = FilterSpec {
            max_info: vec![MaxInfoFilter {
                field: "AF".to_string(),
                max: 0.1,
            }],
            ..empty_spec()
        };
        let sample_index = HashMap::new();
        assert!(spec.keep_record(&record(
            &["chr1", "1", ".", "A", "C", ".", "PASS", "AF=0.01,0.1"],
            &sample_index
        ))?);
        assert!(!spec.keep_record(&record(
            &["chr1", "1", ".", "A", "C", ".", "PASS", "AF=0.01,0.2"],
            &sample_index
        ))?);
        assert!(!spec.keep_record(&record(
            &["chr1", "1", ".", "A", "C", ".", "PASS", "DP=10"],
            &sample_index
        ))?);
        Ok(())
    }

    #[test]
    fn filter_spec_applies_numeric_string_and_missing_expressions() -> Result<()> {
        let spec = FilterSpec {
            expressions: vec![
                parse_expr("gnomAD_AF < 0.01 && DP >= 20")?,
                parse_expr("CLNSIG == 'Pathogenic'")?,
                parse_expr("missing(COMMON)")?,
            ],
            ..empty_spec()
        };
        let sample_index = HashMap::new();
        assert!(spec.keep_record(&record(
            &[
                "chr1",
                "1",
                ".",
                "A",
                "C",
                ".",
                "PASS",
                "gnomAD_AF=0.005;DP=25;CLNSIG=Pathogenic",
            ],
            &sample_index
        ))?);
        assert!(!spec.keep_record(&record(
            &[
                "chr1",
                "1",
                ".",
                "A",
                "C",
                ".",
                "PASS",
                "gnomAD_AF=0.02;DP=25;CLNSIG=Pathogenic",
            ],
            &sample_index
        ))?);
        assert!(!spec.keep_record(&record(
            &[
                "chr1",
                "1",
                ".",
                "A",
                "C",
                ".",
                "PASS",
                "gnomAD_AF=0.005;DP=25;CLNSIG=Pathogenic;COMMON",
            ],
            &sample_index
        ))?);
        Ok(())
    }

    #[test]
    fn filter_spec_applies_boolean_expression_operators() -> Result<()> {
        let spec = FilterSpec {
            expressions: vec![parse_expr(
                "(SOMATIC || CLNSIG == 'Pathogenic') && !COMMON",
            )?],
            ..empty_spec()
        };
        let sample_index = HashMap::new();
        assert!(spec.keep_record(&record(
            &["chr1", "1", ".", "A", "C", ".", "PASS", "SOMATIC;DP=10",],
            &sample_index
        ))?);
        assert!(spec.keep_record(&record(
            &["chr1", "1", ".", "A", "C", ".", "PASS", "CLNSIG=Pathogenic",],
            &sample_index
        ))?);
        assert!(!spec.keep_record(&record(
            &[
                "chr1",
                "1",
                ".",
                "A",
                "C",
                ".",
                "PASS",
                "CLNSIG=Pathogenic;COMMON",
            ],
            &sample_index
        ))?);
        assert!(!spec.keep_record(&record(
            &["chr1", "1", ".", "A", "C", ".", "PASS", "DP=10"],
            &sample_index
        ))?);
        Ok(())
    }

    #[test]
    fn filter_spec_applies_sample_gt_and_dp_predicates() -> Result<()> {
        let spec = FilterSpec {
            sample_has_alt: vec!["Tumor".to_string()],
            sample_gt: vec![SampleGtFilter {
                sample: "Normal".to_string(),
                gt: "0/0".to_string(),
            }],
            sample_min_dp: vec![SampleMinDpFilter {
                sample: "Tumor".to_string(),
                min_dp: 10,
            }],
            ..empty_spec()
        };
        let sample_index = sample_index(&["Tumor", "Normal"]);
        assert!(spec.keep_record(&record(
            &[
                "chr1", "1", ".", "A", "C", ".", "PASS", ".", "GT:DP", "0/1:12", "0/0:20",
            ],
            &sample_index
        ))?);
        assert!(!spec.keep_record(&record(
            &[
                "chr1", "1", ".", "A", "C", ".", "PASS", ".", "GT:DP", "0/0:12", "0/0:20",
            ],
            &sample_index
        ))?);
        assert!(!spec.keep_record(&record(
            &[
                "chr1", "1", ".", "A", "C", ".", "PASS", ".", "GT:DP", "0/1:8", "0/0:20",
            ],
            &sample_index
        ))?);
        Ok(())
    }

    #[test]
    fn filter_spec_applies_group_predicates() -> Result<()> {
        let spec = FilterSpec {
            sample_groups: [
                (
                    "affected".to_string(),
                    vec!["TumorA".to_string(), "TumorB".to_string()],
                ),
                ("controls".to_string(), vec!["Normal".to_string()]),
            ]
            .into_iter()
            .collect(),
            group_any_has_alt: vec!["affected".to_string()],
            group_any_gt: vec![GroupGtFilter {
                group: "affected".to_string(),
                gt: "0/1".to_string(),
            }],
            group_all_min_dp: vec![GroupMinDpFilter {
                group: "controls".to_string(),
                min_dp: 15,
            }],
            ..empty_spec()
        };
        let sample_index = sample_index(&["TumorA", "TumorB", "Normal"]);
        assert!(spec.keep_record(&record(
            &[
                "chr1", "1", ".", "A", "C", ".", "PASS", ".", "GT:DP", "0/0:12", "0/1:14",
                "0/0:20",
            ],
            &sample_index
        ))?);
        assert!(!spec.keep_record(&record(
            &[
                "chr1", "1", ".", "A", "C", ".", "PASS", ".", "GT:DP", "0/0:12", "0/0:14",
                "0/0:20",
            ],
            &sample_index
        ))?);
        assert!(!spec.keep_record(&record(
            &[
                "chr1", "1", ".", "A", "C", ".", "PASS", ".", "GT:DP", "0/0:12", "0/1:14",
                "0/0:10",
            ],
            &sample_index
        ))?);
        Ok(())
    }

    #[test]
    fn filter_vcf_writes_filtered_bgzip_and_index() -> Result<()> {
        let dir = tempdir()?;
        let input = dir.path().join("input.vcf");
        let output = dir.path().join("out.vcf.gz");
        std::fs::write(
            &input,
            "##fileformat=VCFv4.3\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr1\t10\t.\tA\tC\t.\tPASS\tAF=0.05\nchr1\t11\t.\tA\tG\t.\tPASS\tAF=0.2\n",
        )?;
        let spec = FilterSpec {
            max_info: vec![MaxInfoFilter {
                field: "AF".to_string(),
                max: 0.1,
            }],
            ..empty_spec()
        };

        let metrics = filter_vcf(&input, &output, IndexType::Csi, &spec)?;
        assert_eq!(metrics.input_records, 2);
        assert_eq!(metrics.output_records, 1);
        assert!(dir.path().join("out.vcf.gz.csi").exists());

        let mut reader = bgzf::io::Reader::new(File::open(output)?);
        let mut text = String::new();
        reader.read_to_string(&mut text)?;
        assert!(text.contains("chr1\t10\t.\tA\tC\t.\tPASS\tAF=0.05"));
        assert!(!text.contains("chr1\t11\t.\tA\tG"));
        Ok(())
    }

    #[test]
    fn filter_vcf_applies_expression_predicates() -> Result<()> {
        let dir = tempdir()?;
        let input = dir.path().join("input.vcf");
        let output = dir.path().join("out.vcf.gz");
        std::fs::write(
            &input,
            "##fileformat=VCFv4.3\n\
             #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
             chr1\t10\t.\tA\tC\t.\tPASS\tgnomAD_AF=0.005;DP=25;CLNSIG=Pathogenic\n\
             chr1\t11\t.\tA\tG\t.\tPASS\tgnomAD_AF=0.02;DP=25;CLNSIG=Pathogenic\n\
             chr1\t12\t.\tA\tT\t.\tPASS\tgnomAD_AF=0.005;DP=5;CLNSIG=Benign\n",
        )?;
        let spec = FilterSpec {
            expressions: vec![parse_expr(
                "gnomAD_AF < 0.01 && DP >= 20 && CLNSIG != 'Benign'",
            )?],
            ..empty_spec()
        };

        let metrics = filter_vcf(&input, &output, IndexType::Csi, &spec)?;
        assert_eq!(metrics.input_records, 3);
        assert_eq!(metrics.output_records, 1);

        let mut reader = bgzf::io::Reader::new(File::open(output)?);
        let mut text = String::new();
        reader.read_to_string(&mut text)?;
        assert!(text.contains("chr1\t10\t.\tA\tC"));
        assert!(!text.contains("chr1\t11\t.\tA\tG"));
        assert!(!text.contains("chr1\t12\t.\tA\tT"));
        Ok(())
    }

    #[test]
    fn filter_vcf_applies_sample_aware_predicates() -> Result<()> {
        let dir = tempdir()?;
        let input = dir.path().join("input.vcf");
        let output = dir.path().join("out.vcf.gz");
        std::fs::write(
            &input,
            "##fileformat=VCFv4.3\n\
             #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tTumor\tNormal\n\
             chr1\t10\t.\tA\tC\t.\tPASS\tAF=0.05\tGT:DP\t0/1:12\t0/0:20\n\
             chr1\t11\t.\tA\tG\t.\tPASS\tAF=0.05\tGT:DP\t0/0:12\t0/0:20\n\
             chr1\t12\t.\tA\tT\t.\tPASS\tAF=0.05\tGT:DP\t0/1:8\t0/0:20\n",
        )?;
        let spec = FilterSpec {
            sample_has_alt: vec!["Tumor".to_string()],
            sample_min_dp: vec![SampleMinDpFilter {
                sample: "Tumor".to_string(),
                min_dp: 10,
            }],
            sample_gt: vec![SampleGtFilter {
                sample: "Normal".to_string(),
                gt: "0/0".to_string(),
            }],
            ..empty_spec()
        };

        let metrics = filter_vcf(&input, &output, IndexType::Csi, &spec)?;
        assert_eq!(metrics.input_records, 3);
        assert_eq!(metrics.output_records, 1);

        let mut reader = bgzf::io::Reader::new(File::open(output)?);
        let mut text = String::new();
        reader.read_to_string(&mut text)?;
        assert!(text.contains("chr1\t10\t.\tA\tC"));
        assert!(!text.contains("chr1\t11\t.\tA\tG"));
        assert!(!text.contains("chr1\t12\t.\tA\tT"));
        Ok(())
    }
}
