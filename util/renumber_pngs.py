import os

num_digits = 6

lsp = os.popen("ls 0*.png")

i = 0
for f in lsp:
    f = f.rstrip()
    if f.endswith(".png"):
        i_str = str(i)
        name = "renum_" + ("0" * (num_digits - len(i_str))) + i_str + ".png"
        print "mv %s %s" % (f, name)
        os.system("mv %s %s" % (f, name))
        i += 1
