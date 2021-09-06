from Modules.heatmap import heatmap_image, form_examples


xlist = []
ylist = []
modelist = []
RGB_set = [
    [1, 0, 0],
    [0.1, 0.1, 0.1],
    [0.21, 0.75, 0.4],
    [0.7, 0.7, 0.7],
    [0, 0, 1],
]
patterns = ["/", None, "\\", ".", "x"]

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


heatmap_image(dx, dy, xlist, ylist, modelist, RGB_set, patterns)
form_examples(RGB_set, patterns)
print("heat map is formed")
input()
