//! Float File Parser

use core::pbrt::Float;
use pest::iterators::Pair;
use pest::Parser;
use std::fs;
use std::result::Result;

#[derive(Parser)]
#[grammar = "parser/grammar.pest"]
struct FloatParser;

/// Reads a file containing floating point values. If a line contains `#`
/// indicating a comment, the remainder of that line is ignored. Values are
/// read in left to right, top to bottom into a single list. If there are
/// non-comments or numeric values an error is returned.
///
/// NOTE: This will not be very efficient for extremely large files.
///
/// * `path` - Path to file.
pub fn parse_float_file(path: &str) -> Result<Vec<Float>, String> {
    // Load the file and parse the `file` rule.
    let unparsed_file = file_to_string(path)?;
    let file = parse_file_rule(&unparsed_file)?;

    let mut v: Vec<Float> = vec![];
    let mut line_no = 1;

    for nums in file.into_inner() {
        // Parse the `nums` rule for each floating point number.
        match nums.as_rule() {
            Rule::nums => {
                for num in nums.into_inner() {
                    let s = num.as_str();
                    if let Ok(n) = s.parse::<Float>() {
                        v.push(n);
                    } else {
                        return Err(format!(
                            "Error parsing floating point number '{}', line {}.",
                            path, line_no
                        ));
                    }
                }
            }
            Rule::EOI => (), // End of input.
            _ => unreachable!(),
        }

        line_no += 1;
    }

    Ok(v)
}

/// Read the entire file and return its contents as a String.
///
/// * `path` - Path to file.
fn file_to_string(path: &str) -> Result<String, String> {
    match fs::read_to_string(path) {
        Ok(s) => Ok(s),
        _ => Err(format!("Error reading file '{}'", path)),
    }
}

/// Parse the initial `file` rule of the grammar and return the resulting token
/// pairs for remaining rules.
///
/// * `path` - Path to file.
fn parse_file_rule(unparsed_file: &str) -> Result<Pair<'_, Rule>, String> {
    match FloatParser::parse(Rule::file, &unparsed_file) {
        Ok(mut pairs) => Ok(pairs.next().unwrap()), // unwrap `file` rule never fails.
        Err(err) => Err(format!("Error parsing file rule. {}", err)),
    }
}
