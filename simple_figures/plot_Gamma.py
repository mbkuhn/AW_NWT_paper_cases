import numpy as np
import matplotlib.pyplot as plt

N = 100

xlo = 0.
Lg = 4.

def xtilde_g(x):
    xg = 1. - (x - xlo) / Lg
    return xg

def xtilde_r(x, xhi, Lr, Br):
    xr = (x - (xhi - Lr)) / (Lr * Br)
    return xr

def Gamma(xt):
    return 1 - (np.exp(xt**3.5) - 1) / (np.exp(1) - 1)


def main():
    plt.rcParams['text.usetex'] = True

    xg = np.linspace(xlo, xlo + Lg, N)
    Gamma_g = Gamma(xtilde_g(xg))

    xhi = 24.
    Lr = 8.
    Br = 1.
    xr = np.linspace(xhi - Lr, xhi, N)
    Gamma_r = Gamma(xtilde_r(xr, xhi, Lr, Br))
    xr2 = xr
    Br = 2.
    Gamma_r2 = Gamma(xtilde_r(xr2, xhi, Lr, Br))
    xhi = 20.
    Lr = 4.
    xr3 = np.linspace(xhi - Lr, xhi, N)
    Gamma_r3 = Gamma(xtilde_r(xr3, xhi, Lr, Br))

    plt.figure(figsize=(4, 3))
    plt.plot(xg, Gamma_g, "k", label=r'$x_\textrm{lo} = 0, L_g = 4$', linewidth=2)
    plt.xlabel(r'$x$ [m]' )
    plt.ylabel(r'$\Gamma_g$')
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig("plots/Gamma_g.png", dpi=300)

    plt.figure(figsize=(4, 3))
    plt.plot(xr, Gamma_r, "k-",label=r'$x_\textrm{hi} = 24, L_r = 8, B_r = 1$', linewidth=2)
    plt.plot(xr2, Gamma_r2, "b:", label=r'$x_\textrm{hi} = 24, L_r = 8, B_r = 2$', linewidth=2)
    plt.plot(xr3, Gamma_r3, "r--", label=r'$x_\textrm{hi} = 20, L_r = 4, B_r = 2$', linewidth=2)
    plt.xlabel(r'$x$ [m]')
    plt.ylabel(r'$\Gamma_r$')
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig("plots/Gamma_r.png", dpi=300)


if __name__ == "__main__":
    main()
