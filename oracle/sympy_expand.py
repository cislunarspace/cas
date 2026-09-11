#!/usr/bin/env python3
"""sympy oracle：按行读 plain 表达式与有理点，输出 expand 后的精确值。

协议（制表符分隔；与 cas/xtask/src/oracle.rs 对应）：
  请求行: <E|C|F>\t<i>\t<plain 表达式>\t<x y z w 的有理点，形如 3/7 -2/5 ...>
          E = expand 对拍，C = cancel 对拍，F = factor 对拍（点列忽略）
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
        if len(parts) != 4 or parts[0] not in ("E", "C", "F", "D", "T", "S", "A"):
            continue
        tag, i, text, pt = parts[0], parts[1], parts[2], parts[3]
        try:
            # parse_expr 默认不含 convert_xor（^ 会被当作异或），显式替换
            expr = parse_expr(text.replace("^", "**"))
            if tag == "A":
                # 假设对拍：x 声明 positive，simplify 后点值（正点）
                from sympy import Symbol, simplify as sym_simplify
                xp = Symbol("x", positive=True)
                xnum, xden = pt.split()[0].split("/")
                env = {xp: Rational(int(xnum), int(xden))}
                expr_a = expr.subs(VARS[0], xp)
                val = sym_simplify(expr_a).subs(env)
                fv = val.evalf(30)
                if not fv.is_number or not fv.is_finite:
                    print(f"V\t{i}\tERR", flush=True)
                    continue
                print(f"V\t{i}\tF:{fv}", flush=True)
                continue
            if tag in ("D", "T", "S"):
                env = {}
                for s_, p_ in zip(VARS, pt.split()):
                    n_, d_ = p_.split("/")
                    env[s_] = Rational(int(n_), int(d_))
                # diff / taylor / simplify 对拍：对表达式做对应变换后在点求值
                from sympy import diff as sym_diff
                from sympy.series import series as sym_series
                xv = VARS[0]  # 主变元固定 x
                if tag == "D":
                    val = sym_diff(expr, xv).subs(env)
                elif tag == "T":
                    # x=0 处 5 阶（与 xtask 端约定一致）
                    ser = sym_series(expr, xv, 0, 6).removeO()
                    val = ser.subs(env)
                else:
                    from sympy import simplify as sym_simplify
                    val = sym_simplify(expr).subs(env)
                # 三角在整数点的值是符号的：浮点对拍（30 位精度，容差在
                # xtask 端控制）
                fv = val.evalf(30)
                if not fv.is_number or not fv.is_finite:
                    print(f"V\t{i}\tERR", flush=True)
                    continue
                print(f"V\t{i}\tF:{fv}", flush=True)
                continue
            if tag == "F":
                # 因式分解对拍：返回 重数:次数 多重集（排序后序列化）；不需点
                from sympy import factor_list, total_degree
                _, facs = factor_list(expr)
                ms = sorted(f"{m}:{total_degree(f)}" for f, m in facs)
                print(f"V\t{i}\t{','.join(ms)}", flush=True)
                continue
            env = {}
            for s, p in zip(VARS, pt.split()):
                n, d = p.split("/")
                env[s] = Rational(int(n), int(d))
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
