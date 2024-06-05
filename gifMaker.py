import imageio


directory = "VelocityProfile_30000_it"
maxIter = 30000
plots = 100

frames = []
t = maxIter//plots

print("Making Gif ...")
for i in range(t):
    # print(i)
    image = imageio.v2.imread(f"./" + directory + "/profile_" + str(i) + ".png")
    frames.append(image)

imageio.mimsave("./" + directory + "/_" + directory + ".gif", frames, duration = 80)

print("\nDone.")