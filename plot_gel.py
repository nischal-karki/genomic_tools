from matplotlib import pyplot as plt
from gen_gel_map import gen_frag_sizes

def plot_gel(ax,file_path, restriction_enzymes, index=0):
    gfs = gen_frag_sizes(file_path, restriction_enzymes)
    gfs.gen_frag_sizes()
    sizes = gfs.frag_sizes
    plotted = ax.scatter([index]*len(sizes), sizes, s=100, marker="_")
    ax.set_yscale("log")
    for s in sizes:
        ax.text(index, s, str(s), ha="center", va="bottom")
    return plotted
if __name__ == "__main__":
    import sys
    import glob
    import os

    if "--label" in sys.argv or "-l" in sys.argv:
        label = sys.argv.index("--label") if "--label" in sys.argv else sys.argv.index("-l")
        sys.argv.pop(label)
        label = sys.argv.pop(label)
    else:
        label = None
    files = None
    if "--start_enzyme_here" in sys.argv:
        start = sys.argv.index("--start_enzyme_here")
        sys.argv.pop(start)
    elif "-e" in sys.argv:
        start = sys.argv.index("-e")
        sys.argv.pop(start)
    else:
        files = [sys.argv.pop(i) for i,arg in enumerate(sys.argv) if arg.endswith((".fasta",".fa",".gb",".gbk"))]
        start = 1
    files = sys.argv[1:start] if files is None else files
    enzymes = sys.argv[start:]

    if "--glob" == sys.argv[1] or "-g" == sys.argv[1]:
        sys.argv.pop(1)
        files = glob.glob(sys.argv[1])
        enzymes = sys.argv[2:]

    fig = plt.figure(dpi=300)
    axes = fig.add_subplot(111)
    for i,file in enumerate(files):
        if label is None:
            label = os.path.basename(file)
        else:
            count = 1
            while label in os.listdir():
                label = f"{label}_{count}"
                count += 1
        plot_gel(axes, file, enzymes,index=i).set_label(label)
    plt.legend(loc="upper center",ncol=3, bbox_to_anchor=(0.5,1.2),fontsize=6)
    [spine.set_visible(False) for spine in axes.spines.values()]
    axes.set_xticks([])
    axes.set_yticks([])
    axes.minorticks_off()
    plt.tight_layout()
    plt.savefig("gel_map.png")
    plt.show()