local: generate_figures externalize compile_paper
push: generate_figures externalize compile_paper push_paper

generate_figures:
	julia scripts_jl/paper_figures/make_figures.jl

externalize:
	cd externalization; python externalize.py

compile_paper:
	cd paper; latexmk -pdf main.tex

push_paper:
	cd paper; git add *; git commit -m "make"; git push origin master