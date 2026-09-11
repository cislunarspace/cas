#!/usr/bin/env python3
"""sympy oracle：按行读 plain 表达式与有理点，输出 expand 后的精确值。

协议（制表符分隔；与 cas/xtask/src/oracle.rs 对应）：
  请求行: <E|C>\t<i>\t<plain 表达式>\t<x y z w 的有理点，形如 3/7 -2/5 ...>
          E = expand 对拍，C = cancel 对拍
  响应行: V\t<i>\t<num>/<den>；非有限值（zoo/oo/nan）输出 ERR

依赖：sympy（版本锁定见 requirements.txt）。
"""
import sys
from sympy import Rational, cancel, expand, symbols
from sympy.parsing.sympy_parser import parse_expr

VARS = symbols("x y z w")


def main() -> None:
    for line in sys.stdin:
        parts = line.rstrip("\n").split("\t")
        if len(parts) != 4 or parts[0] not in ("E", "C"):
            continue
        tag, i, text, pt = parts[0], parts[1], parts[2], parts[3]
        try:
            env = {}
            for s, p in zip(VARS, pt.split()):
                n, d = p.split("/")
                env[s] = Rational(int(n), int(d))
            # parse_expr 默认不含 convert_xor（^ 会被当作异或），显式替换
            expr = parse_expr(text.replace("^", "**"))
            val = (cancel(expr) if tag == "C" else expand(expr)).subs(env)
            if not val.is_Rational:
                print(f"V\t{i}\tERR", flush=True)
                continue
            n, d = val.as_numer_denom()
            print(f"V\t{i}\t{n}/{d}", flush=True)
        except Exception:
            print(f"V\t{i}\tERR", flush=True)


if __name__ == "__main__":
    main()
