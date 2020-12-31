//! Common

use super::Float;
use nom::branch::alt;
use nom::bytes::complete::take_till;
use nom::character::complete::{char, one_of, space1};
use nom::combinator::{opt, recognize};
use nom::multi::{many0, many1, separated_list0};
use nom::sequence::{preceded, terminated, tuple};
use nom::IResult;
use std::fs::File;
use std::io::{BufRead, BufReader};

/// Reads a file containing floating point values. If a line contains `#`
/// indicating a comment, the remainder of that line is ignored. Values are
/// read in left to right, top to bottom into a single list. If there are
/// non-comments or numeric values an error is returned.
pub fn read_float_file(path: &str) -> Result<Vec<Float>, String> {
    match File::open(path) {
        Ok(file) => {
            let reader = BufReader::new(file);
            let mut v: Vec<Float> = vec![];

            for (index, line) in reader.lines().enumerate() {
                match line {
                    Ok(s) => {
                        let d = process_line(index, &s)?;
                        v.extend(d.iter());
                    }
                    Err(err) => {
                        return Err(format!("Error reading {}. {}.", path, err));
                    }
                }
            }

            Ok(v)
        }
        Err(err) => Err(format!("Error reading {}. {}.", path, err)),
    }
}

/// Process a single line. If the line cannot be parsed entirely an error is
/// returned.
///
/// * `index` - The line number (starting at 0).
/// * `line`  - The contents of the line.
fn process_line(index: usize, line: &str) -> Result<Vec<Float>, String> {
    match parse_line(line) {
        Ok((remaining, data)) => {
            if remaining.len() == 0 {
                let mut v: Vec<Float> = vec![];
                for s in data {
                    match s.parse::<Float>() {
                        Ok(d) => v.push(d),
                        Err(err) => {
                            return Err(format!("Error parsing Float value {}. {}.", s, err))
                        }
                    }
                }
                Ok(v)
            } else {
                Err(format!(
                    "Could not parse remainder of line {}, `{}`.",
                    index + 1,
                    remaining
                ))
            }
        }
        Err(err) => Err(format!("Error parsing line {}. {}.", index + 1, err)),
    }
}

/// Strips any comments and parses floating point values (including decimals) and
/// returns a list of parsed values and remaining unparsed string.
///
/// * `input` - The input string.
fn parse_line(input: &str) -> IResult<&str, Vec<&str>> {
    let (_comment, leading) = take_till_comment(input)?;
    floats(leading.trim()) // trim needed because separated_list0 doesn't like leading whitepsace.
}

/// Parses the string upto a `#` symbol and returns a tuple containing the
/// (comment, leading string).
///
/// * `input` - The input string.
fn take_till_comment(input: &str) -> IResult<&str, &str> {
    take_till(|c| c == '#')(input)
}

/// Parses a sequence of floating point values (including decimal numbers)
/// and returns a list of parsed values and remaining unparsed string.
///
/// * `input` - The input string.
fn floats(input: &str) -> IResult<&str, Vec<&str>> {
    separated_list0(space1, alt((float, decimal)))(input)
}

/// Parses a decimal number and produces the remainder of the unparsed string.
/// If no decimal number is found, the remainder is the whole input string.
///
/// * `input` - The input string.
fn decimal(input: &str) -> IResult<&str, &str> {
    recognize(many1(terminated(one_of("0123456789"), many0(char('_')))))(input)
}

/// Parses a floating point number and produces the remainder of the unparsed
/// string. If no floating point number is found, the remainder is the whole input
/// string.
///
/// * `input` - The input string.
fn float(input: &str) -> IResult<&str, &str> {
    alt((
        // Case one: .42
        recognize(tuple((
            char('.'),
            decimal,
            opt(tuple((one_of("eE"), opt(one_of("+-")), decimal))),
        ))),
        // Case two: 42e42 and 42.42e42
        recognize(tuple((
            decimal,
            opt(preceded(char('.'), decimal)),
            one_of("eE"),
            opt(one_of("+-")),
            decimal,
        ))),
        // Case three: 42. and 42.42
        recognize(tuple((decimal, char('.'), opt(decimal)))),
    ))(input)
}
