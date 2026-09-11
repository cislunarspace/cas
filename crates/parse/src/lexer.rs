//! 手写词法：ASCII 标识符、十进制整数、浮点（必须含小数点，可带指数）。
//!
//! 注意：`1e3` 不是浮点——指数记号仅随小数点出现；裸 `e3` 按标识符处理，
//! 因此 `1e3` 会解析失败（期望运算符）。这是刻意的文法决定，保证词法无歧义。

use cas_domain::Integer;

#[derive(Clone, Debug)]
pub(crate) enum Tok<'s> {
    Int(Integer),
    Float(f64),
    Ident(&'s str),
    Plus,
    Minus,
    Star,
    Slash,
    Caret,
    LParen,
    RParen,
    Comma,
}

#[derive(Clone, Debug)]
pub(crate) struct Token<'s> {
    pub tok: Tok<'s>,
    pub pos: usize,
}

pub(crate) fn lex(src: &str) -> Result<Vec<Token<'_>>, super::ParseError> {
    let b = src.as_bytes();
    let mut i = 0;
    let mut out = Vec::new();
    while i < b.len() {
        let c = b[i];
        if matches!(c, b' ' | b'\t' | b'\r' | b'\n') {
            i += 1;
            continue;
        }
        let start = i;
        if c.is_ascii_digit() {
            while i < b.len() && b[i].is_ascii_digit() {
                i += 1;
            }
            if i < b.len() && b[i] == b'.' {
                i += 1;
                let frac_start = i;
                while i < b.len() && b[i].is_ascii_digit() {
                    i += 1;
                }
                if i == frac_start {
                    return Err(super::ParseError::new(frac_start, "小数点后缺数字"));
                }
                if i < b.len() && (b[i] == b'e' || b[i] == b'E') {
                    let e_pos = i;
                    i += 1;
                    if i < b.len() && (b[i] == b'+' || b[i] == b'-') {
                        i += 1;
                    }
                    let exp_start = i;
                    while i < b.len() && b[i].is_ascii_digit() {
                        i += 1;
                    }
                    if i == exp_start {
                        return Err(super::ParseError::new(e_pos, "指数缺数字"));
                    }
                }
                let v: f64 = src[start..i]
                    .parse()
                    .map_err(|_| super::ParseError::new(start, "浮点字面量非法"))?;
                out.push(Token {
                    tok: Tok::Float(v),
                    pos: start,
                });
            } else {
                let v = Integer::parse(&src[start..i])
                    .ok_or_else(|| super::ParseError::new(start, "整数字面量解析失败"))?;
                out.push(Token {
                    tok: Tok::Int(v),
                    pos: start,
                });
            }
            continue;
        }
        if c.is_ascii_alphabetic() || c == b'_' {
            i += 1;
            while i < b.len() && (b[i].is_ascii_alphanumeric() || b[i] == b'_') {
                i += 1;
            }
            out.push(Token {
                tok: Tok::Ident(&src[start..i]),
                pos: start,
            });
            continue;
        }
        let tok = match c {
            b'+' => Tok::Plus,
            b'-' => Tok::Minus,
            b'*' => Tok::Star,
            b'/' => Tok::Slash,
            b'^' => Tok::Caret,
            b'(' => Tok::LParen,
            b')' => Tok::RParen,
            b',' => Tok::Comma,
            _ => {
                let ch = src[start..].chars().next().unwrap_or('?');
                return Err(super::ParseError::new(start, &format!("非法字符 {ch:?}")));
            }
        };
        i += 1;
        out.push(Token { tok, pos: start });
    }
    Ok(out)
}
