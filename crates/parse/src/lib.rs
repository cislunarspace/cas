//! 解析层：plain 文法（与 cas-print 的 plain 输出互为往返，D6 承诺）。
//!
//! 文法（优先级 `^` > 一元 ± > `*` `/` > `+` `-`；`^` 右结合，指数可带符号）：
//!
//! ```text
//! expr    := term (('+' | '-') term)*
//! term    := unary (('*' | '/') unary)*
//! unary   := ('+' | '-') unary | postfix
//! postfix := atom ('^' unary)?
//! atom    := int | float | ident | ident '(' expr (',' expr)* ')' | '(' expr ')'
//! ```
//!
//! 标识符后随 `(` 视为函数应用（参数保持书写次序），否则为符号。
//! 构造经 `Context` 的规范形入口，解析结果即刻规范化。

mod lexer;

use cas_expr::{Context, Expr};
use lexer::{Tok, Token};
use std::error::Error;
use std::fmt;

/// 解析错误：`pos` 为字节偏移。
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ParseError {
    pub pos: usize,
    pub msg: String,
}

impl ParseError {
    fn new(pos: usize, msg: &str) -> Self {
        ParseError {
            pos,
            msg: msg.to_string(),
        }
    }
}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "位置 {}: {}", self.pos, self.msg)
    }
}

impl Error for ParseError {}

pub fn parse(ctx: &Context, src: &str) -> Result<Expr, ParseError> {
    let toks = lexer::lex(src)?;
    let mut p = Parser {
        ctx,
        toks,
        idx: 0,
        src_len: src.len(),
    };
    let e = p.expr()?;
    match p.peek() {
        Some(t) => Err(ParseError::new(t.pos, "存在多余记号")),
        None => Ok(e),
    }
}

struct Parser<'a> {
    ctx: &'a Context,
    toks: Vec<Token<'a>>,
    idx: usize,
    src_len: usize,
}

impl<'a> Parser<'a> {
    fn peek(&self) -> Option<Token<'a>> {
        self.toks.get(self.idx).cloned()
    }

    fn expr(&mut self) -> Result<Expr, ParseError> {
        let mut l = self.term()?;
        loop {
            let plus = match self.peek().map(|t| t.tok) {
                Some(Tok::Plus) => true,
                Some(Tok::Minus) => false,
                _ => break,
            };
            self.idx += 1;
            let r = self.term()?;
            l = if plus {
                self.ctx.add(&[l, r])
            } else {
                let neg = self.ctx.mul(&[r, self.ctx.int(-1)]);
                self.ctx.add(&[l, neg])
            };
        }
        Ok(l)
    }

    fn term(&mut self) -> Result<Expr, ParseError> {
        let mut l = self.unary()?;
        loop {
            let mul = match self.peek().map(|t| t.tok) {
                Some(Tok::Star) => true,
                Some(Tok::Slash) => false,
                _ => break,
            };
            self.idx += 1;
            let r = self.unary()?;
            l = if mul {
                self.ctx.mul(&[l, r])
            } else {
                let inv = self.ctx.pow(&r, &self.ctx.int(-1));
                self.ctx.mul(&[l, inv])
            };
        }
        Ok(l)
    }

    fn unary(&mut self) -> Result<Expr, ParseError> {
        match self.peek().map(|t| t.tok) {
            Some(Tok::Minus) => {
                self.idx += 1;
                let e = self.unary()?;
                Ok(self.ctx.mul(&[e, self.ctx.int(-1)]))
            }
            Some(Tok::Plus) => {
                self.idx += 1;
                self.unary()
            }
            _ => self.postfix(),
        }
    }

    fn postfix(&mut self) -> Result<Expr, ParseError> {
        let a = self.atom()?;
        if matches!(self.peek().map(|t| t.tok), Some(Tok::Caret)) {
            self.idx += 1;
            let e = self.unary()?;
            Ok(self.ctx.pow(&a, &e))
        } else {
            Ok(a)
        }
    }

    fn atom(&mut self) -> Result<Expr, ParseError> {
        let Some(t) = self.peek() else {
            return Err(ParseError::new(self.src_len, "期望表达式，遇到结尾"));
        };
        match t.tok {
            Tok::Int(ref v) => {
                self.idx += 1;
                Ok(self.ctx.integer(v))
            }
            Tok::Float(v) => {
                self.idx += 1;
                Ok(self.ctx.float(v))
            }
            Tok::Ident(name) => {
                self.idx += 1;
                if matches!(self.peek().map(|t| t.tok), Some(Tok::LParen)) {
                    self.idx += 1;
                    let mut args = vec![self.expr()?];
                    while matches!(self.peek().map(|t| t.tok), Some(Tok::Comma)) {
                        self.idx += 1;
                        args.push(self.expr()?);
                    }
                    let close = self.peek();
                    match close.as_ref().map(|t| &t.tok) {
                        Some(Tok::RParen) => {
                            self.idx += 1;
                            Ok(self.ctx.call(name, &args))
                        }
                        _ => Err(ParseError::new(
                            close.as_ref().map(|t| t.pos).unwrap_or(self.src_len),
                            "期望 ')' 结束函数参数",
                        )),
                    }
                } else {
                    Ok(self.ctx.sym(name))
                }
            }
            Tok::LParen => {
                self.idx += 1;
                let e = self.expr()?;
                match self.peek().map(|t| t.tok) {
                    Some(Tok::RParen) => {
                        self.idx += 1;
                        Ok(e)
                    }
                    _ => Err(ParseError::new(self.idx_at(), "期望 ')'")),
                }
            }
            _ => Err(ParseError::new(t.pos, "期望表达式")),
        }
    }

    fn idx_at(&self) -> usize {
        self.toks
            .get(self.idx)
            .map(|t| t.pos)
            .unwrap_or(self.src_len)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use cas_print::plain;

    fn round(src: &str, expect: &str) {
        let ctx = Context::new();
        let e = parse(&ctx, src).unwrap_or_else(|err| panic!("{src} 解析失败: {err}"));
        assert_eq!(plain(&ctx, &e), expect, "输入: {src}");
    }

    #[test]
    fn 优先级与结合性() {
        round("1 + 2 * 3", "7");
        round("x^2*y", "x^2*y");
        round("-x^2", "-x^2");
        round("x^-2", "x^-2");
        round("2^-2", "1/4");
        round("x/y/z", "x*y^-1*z^-1"); // ((x/y)/z) = x · y⁻¹ · z⁻¹
        round("2^3^2", "512"); // 右结合 2^(3^2)，数值幂折叠
        round("(x + y)*x", "x*(x + y)");
        round("7/2", "7/2");
    }

    #[test]
    fn 函数与符号() {
        round("sin(x)", "sin(x)");
        round("sin(x)^2", "sin(x)^2");
        round("sin (x)", "sin(x)");
        round("atan2(y, x)", "atan2(y, x)");
        round("f(f(x))", "f(f(x))");
        round("sin", "sin"); // 无括号：普通符号
        round("b_2 + a", "a + b_2");
    }

    #[test]
    fn 大整数字面量() {
        round(
            "123456789012345678901234567890 + 1",
            "123456789012345678901234567891",
        );
    }

    #[test]
    fn 错误报告() {
        let ctx = Context::new();
        assert!(parse(&ctx, "").is_err());
        assert!(parse(&ctx, "2x").is_err()); // 缺运算符
        assert!(parse(&ctx, "1e3").is_err()); // 指数须随小数点
        assert!(parse(&ctx, "(x").is_err());
        assert!(parse(&ctx, "x)").is_err());
        assert!(parse(&ctx, "f()").is_err()); // 函数至少一个参数
        assert!(parse(&ctx, "1.").is_err());
        assert!(parse(&ctx, "1.5e+").is_err());
        assert!(parse(&ctx, "$").is_err());
        // 位置信息
        let err = parse(&ctx, "x + $").unwrap_err();
        assert_eq!(err.pos, 4);
    }
}
