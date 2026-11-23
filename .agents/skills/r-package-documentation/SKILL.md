---
name: r-package-documentation
description: Write and maintain R package documentation, including Markdown files. The current workflow covers READMEs. Use it to create or revise a README, trace its source, update badges, or migrate equations to native GitHub math.
---

# R package documentation

Apply the `ste-writing` skill to package prose. Use the `r-package-development` skill for package changes and checks.
The workflow below covers READMEs. For other documentation, use the source and writing guidance in `ste-writing`.

## Find the README source

1. Read DESCRIPTION, the README, and the repository instructions.
2. Trace the README build through scripts, configuration, and generated-file notices. A vignette can supply both README.md and an HTML article.
3. Edit the maintained source. For a hand-written README, edit README.md directly.
4. Identify the audience, the package purpose, and the installation method. Check each example against the current R interface.
5. Keep technical terms consistent with the package. Explain the inputs, the result, and what to look for in a plot or table.

## Add badges

Use the repository owner, repository name, workflow filename, and case-sensitive package name in these templates.
Check the workflow under `.github/workflows/`. Use the CRAN badges for a published package.

```markdown
[![R-CMD-check](https://github.com/OWNER/REPO/actions/workflows/WORKFLOW/badge.svg)](https://github.com/OWNER/REPO/actions/workflows/WORKFLOW)
[![CRAN version](https://www.r-pkg.org/badges/version/PACKAGE)](https://CRAN.R-project.org/package=PACKAGE)
[![CRAN downloads](https://cranlogs.r-pkg.org/badges/last-month/PACKAGE)](https://cranlogs.r-pkg.org/)
```

Check the CRAN version image and destination against `usethis::use_cran_badge()`.
The helper can edit the README. Inspect its implementation when the repository generates the README from another source.

## Migrate equation images

1. Recover the original TeX from the maintained source. For a hand-written image link, decode the math query once and compare it with the alt text.
2. Use `$...$` for inline math. Put `$$` on separate lines around display math in the source.
3. For several aligned rows, put `\begin{aligned}` and `\end{aligned}` inside the display. Preserve alignment points, row breaks, subscripts, and grouping.
   If GitHub turns a continuation line beginning with `+ ` into a list item, move the plus sign to the preceding line.
4. Inspect the writer options. Remove legacy `--webtex` arguments. With current rmarkdown and Pandoc, use `rmarkdown::github_document(math_method = "default")`.
5. Regenerate each tracked output with the repository build command. Check the generated Markdown for native delimiters and remaining equation-image URLs.
6. Check the GitHub preview for missing symbols, visible TeX, broken alignment, and clipped rows. Check ordinary figures, badges, and heading links separately.
7. If a reviewer finds a rendering issue, reproduce it in GitHub. Add the verified fix to this skill after checking the revised preview.

[GitHub math syntax](https://docs.github.com/en/get-started/writing-on-github/working-with-advanced-formatting/writing-mathematical-expressions) · [rmarkdown GitHub output](https://pkgs.rstudio.com/rmarkdown/reference/github_document.html)
