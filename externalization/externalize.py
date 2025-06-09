import os
import subprocess
import shutil

# Look for all .tex files generated from pgfplotsx Julia scripts
with open("externalize_subfile.tex", "w") as f:
    for fname in filter(lambda x: x.endswith(".tex"), os.listdir("figures_tex")):
        f.write(r"\includetikz{" + fname[:-4] + "}\n")
        this_pdf = f"figures_pdf/{fname[:-4]}.pdf"
        if os.path.exists(this_pdf):
            os.remove(this_pdf)

# externalize_subfile gets inputted into image_externalization, so we now compile:     
subprocess.run(["pdflatex", "-shell-escape", "-interaction=batchmode", "image_externalization.tex"])

# this creates pdfs in figures_pdf, which we move to the paper subdirectory
for fname in filter(lambda x: x.endswith(".pdf"), os.listdir("figures_pdf")):
    shutil.copy(f"figures_pdf/{fname}", f"../paper/figures/{fname}")