# Ascites

To render R Markdown file
1. Change front matter to:
```
output: html_document
```
2. Run R from command line, where PATH is the path to the R Markdown file:
```
library(rmarkdown)
render(PATH)
```