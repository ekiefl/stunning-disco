import matplotlib.pyplot as plt
import matplotlib

m = matplotlib.cm.get_cmap("tab20")
for i in range(20):
    rgb = m(i)[:3]
    print(matplotlib.colors.rgb2hex(rgb))
