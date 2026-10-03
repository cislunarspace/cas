# 贡献指南

欢迎 Bug 报告、功能建议与代码贡献。技术讨论就事论事;提交即表示同意以与仓库一致的许可授权你的贡献。

## 提 Issue
先搜索避免重复,按模板填写;标题以 [BUG] / [FEAT] 开头,后接一句中文结论句。AI 生成的 issue 与 PR,标题最前面加 [AI Generated]。

## 提 PR
1. 先开 Issue,PR 描述写 `Fixes #NN`;纯文档小修可不挂。
2. 从默认分支切 `fix/<简述>` 或 `feat/<简述>` 分支。
3. 本地验证通过(见各仓库 CONTRIBUTING),改动保持最小,不含无关格式化。
4. CI 通过后 squash 合并。

## 本地验证
各仓库验证命令以其自身 CONTRIBUTING.md 或 CI 工作流为准。

## 本地验证
与 CI 一致:cargo fmt --all -- --check;cargo clippy --workspace --all-targets --locked -- -D warnings;cargo test --workspace --locked。
