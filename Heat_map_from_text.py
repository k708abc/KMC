from Modules.heatmap import heatmap_image


xlist = []
ylist = []
modelist = []

with open("growth_modes.txt", "r") as f:
    for line in f:
        line_s = line.split("\t")
        if line_s[0] == "dx":
            dx = float(line_s[1])
        elif line_s[0] == "dy":
            dy = float(line_s[1])
        else:
            xlist.append(float(line_s[0]))
            ylist.append(float(line_s[1]))
            modelist.append(line_s[2].replace("\n", ""))


heatmap_image(dx, dy, xlist, ylist, modelist)
