# -*- coding: utf-8 -*-
import argparse, json, pylab, sys, math, re, Image, PIL
import numpy
import PIL.Image
import PIL.ImageDraw
import PIL.ImageFont


def draw_text_at_center(img, text):
  draw = PIL.ImageDraw.Draw(img)
  draw.font = PIL.ImageFont.truetype(
    "/usr/share/fonts/truetype/freefont/FreeMono.ttf", 50)

  img_size = numpy.array(img.size)
  txt_size = numpy.array(draw.font.getsize(text))
  pos = (img_size - txt_size) / 2

  draw.text(pos, text, (0, 0, 0))


def make_image(i, screen, bgcolor, filename):
    img = Image.new('RGB', screen, bgcolor)
    draw_text_at_center(img, str(i))
    img.save(filename)

if __name__ == '__main__':
    for i in range(0, 98):
        screen = (1574,70)
        bgcolor=(255,255,255)
        filename = "psi_real_%08d_num.png" % i
        make_image(i + 1700, screen, bgcolor, filename)
