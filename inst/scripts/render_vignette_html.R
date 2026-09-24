# Render vignette using knitr / markdown into standalone HTML
cat("Rendering vignettes/spconform-intro.Rmd...\n")

if (requireNamespace("knitr", quietly = TRUE)) {
  # Knit to markdown first
  md_file <- "vignettes/spconform-intro.md"
  html_file <- "vignettes/spconform-intro.html"
  
  knitr::knit("vignettes/spconform-intro.Rmd", output = md_file)
  
  if (requireNamespace("markdown", quietly = TRUE)) {
    markdown::mark_html(md_file, output = html_file)
  } else {
    # Generate standalone luxury HTML wrapper
    md_content <- readLines(md_file, warn = FALSE)
    # Simple HTML conversion
    html_lines <- c(
      "<!DOCTYPE html>",
      "<html><head><meta charset='utf-8'>",
      "<title>Introduction to spconform</title>",
      "<link rel='stylesheet' href='https://cdn.jsdelivr.net/npm/github-markdown-css@5/github-markdown.min.css'>",
      "<script src='https://cdn.jsdelivr.net/npm/marked/marked.min.js'></script>",
      "<style>",
      "body { box-sizing: border-box; min-width: 200px; max-width: 980px; margin: 0 auto; padding: 45px; background: #0b0f19; color: #f3f4f6; font-family: -apple-system,BlinkMacSystemFont,'Segoe UI',Helvetica,Arial,sans-serif; }",
      ".markdown-body { background: #111827; color: #e5e7eb; padding: 40px; border-radius: 16px; border: 1px solid rgba(255,255,255,0.1); box-shadow: 0 10px 30px rgba(0,0,0,0.5); }",
      ".markdown-body a { color: #60a5fa; }",
      ".markdown-body code { background: #1f2937; color: #38bdf8; }",
      ".markdown-body pre { background: #1f2937; border: 1px solid rgba(255,255,255,0.08); }",
      "</style>",
      "</head><body>",
      "<article id='content' class='markdown-body'></article>",
      "<script>",
      "const rawMd = " , jsonlite::toJSON(paste(md_content, collapse = "\n")), ";",
      "document.getElementById('content').innerHTML = marked.parse(rawMd);",
      "</script>",
      "</body></html>"
    )
    writeLines(html_lines, html_file)
  }
  cat("Vignette rendered successfully to: ", html_file, "\n")
}
