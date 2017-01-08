# -*- coding: utf-8 -*-
# ffmpegの連番画像からの動画生成のために,
# lsで出力される順番で.pngファイルのファイル名を振り直す.
import os

num_digits = 6

lsp = os.popen("ls 0*.png")
i = 0
for f in lsp:
    f = f.rstrip()
    if f.endswith(".png"):
        i_str = str(i)
        name = "a" + ("0" * (num_digits - len(i_str))) + i_str + ".png"
        print "mv %s %s" % (f, name)
        os.system("mv %s %s" % (f, name))
        i += 1
