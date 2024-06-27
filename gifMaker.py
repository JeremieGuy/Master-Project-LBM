import imageio


directory = "Flow_50000_it"
maxIter = 50000
plots = 20

frames = []
t = maxIter//plots

print("Making Gif ...")
for i in range(t):
    # print(i)
    num = "{0:0=5d}".format(i)
    image = imageio.v2.imread(f"./" + directory + "/fluid_" + str(num) + ".png")
    frames.append(image)

imageio.mimsave("./" + directory + "/_" + directory + ".gif", frames, duration = 80)

print("\nDone.")