library("bookdown")
create_gitbook(".")
file.create(".nojekyll")
rmarkdown::render_site(encoding = 'UTF-8')