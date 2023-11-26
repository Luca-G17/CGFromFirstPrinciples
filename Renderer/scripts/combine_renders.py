import os
from os import listdir
from os.path import isfile, join

path = "Renders/" 

def CombineRenders(Renders):
    files = [f for f in listdir(path + Renders[0]) if isfile(join(path + Renders[0], f))]
    n = int(max(files, key=lambda f : int(f.split(".bmp")[0])).split(".bmp")[0])
    files2 = [f for f in listdir(path + Renders[1]) if isfile(join(path + Renders[1], f))]
    for f in files2:
        os.rename(path + Renders[1] + f, f"{path}/{Renders[0]}/{(int(f.split('.bmp')[0]) + n + 1)}.bmp")


def ZeroRenders(Render):
    files = [f for f in listdir(Render) if isfile(join(Render, f))]
    n = int(min(files, key=lambda f : int(f.split(".bmp")[0])).split(".bmp")[0]) - 1
    os.mkdir(Render + "/Zeroed")
    for f in files:
        os.rename(Render + f, f"{Render}/Zeroed/{int(f.split('.bmp')[0]) - n}.bmp")

Renders = ["VideoP1/", "VideoP2/Zeroed/"]
CombineRenders(Renders)

# ZeroRenders("Renders/VideoP2/")