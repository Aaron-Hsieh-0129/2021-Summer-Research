import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import seaborn as sns
import scipy.stats as stats
import time

def getData(file):
	nx, nz = int(150000 / 250), int(15000 / 250) 
	output = np.zeros([nx, nz])
	for k in range(nz):
		for i in range(nx):
			output[i][k] = file.pop(0)

	return output

def getData2(file):
	nx, nz = int(150000 / 250), int(15000 / 250) 
	output = np.zeros([nx, nz])
	for k in range(nz):
		for i in range(nx):
			output[i][k] = file.pop(0)
			if output[i][k] < 0.:
				output[i][k] = 0.;

	return output


dt = 0.1
nx, nz = int(150000 / 250), int(15000 / 250)

count = 1
for t in range(0, 50001, 250):
	
	th = np.loadtxt("outputs/th/th_" + str(t) + ".txt").reshape([nx, nz], order='F')
	th_flip = th.swapaxes(0, 1)[:, 250:350+1]
	plt.figure(figsize=(10, 6))
	plt.tight_layout()
	plt.title(f"t = {t * dt} s", fontsize=14)
	plt.xlabel("x [km]")
	plt.ylabel("z [km]")
	plt.xticks(np.linspace(0, 100, 7), ["60", "65", "70", "75", "80", "85", "90"])
	plt.yticks([0, 10, 20, 30, 40, 50, 59], ["0", "2.5", "5", "7.5", "10", "12.5", "15"])
	CF = plt.contourf(th_flip, levels=np.arange(-16, 16+2,2), extend='both', cmap=cm.twilight_shifted)
	cbar = plt.colorbar(pad=0.05)
	cbar.set_ticks(np.arange(-16, 16+4, 4))
	cbar.set_label(r"$\theta^{'}$ [K]")

	X, Y = np.meshgrid(np.arange(101), np.arange(60))
	u = np.loadtxt("outputs/u/u_" + str(t) + ".txt").reshape([nx, nz], order='F')
	u_flip = u.swapaxes(0, 1)[:, 250:350+1]
	w = np.loadtxt("outputs/w/w_" + str(t) + ".txt").reshape([nx, nz], order='F')
	w_flip = w.swapaxes(0, 1)[:, 250:350+1]
	Q = plt.quiver(X[::2, ::2], Y[::2, ::2], u_flip[::2, ::2], w_flip[::2, ::2], angles='xy', units="width", scale=300)
	qk = plt.quiverkey(Q, 0.7, 0.9, 10, r'$10 \frac{m}{s}$', labelpos='E', coordinates='figure')

	qc = np.loadtxt("outputs/qc/qc_" + str(t) + ".txt").reshape([nx, nz], order='F')
	qc_flip = qc.swapaxes(0, 1)[:, 250:350+1] * 1000
	CS = plt.contour(qc_flip, levels=[0.0001, 1, 2, 3], colors='yellow', linewidths=1.25)
	plt.clabel(CS, inline=1, fontsize=8, fmt="%1.0f")

	qr = np.loadtxt("outputs/qr/qr_" + str(t) + ".txt").reshape([nx, nz], order='F')
	qr_flip = qr.swapaxes(0, 1)[:, 250:350+1]
	x, z = [], []
	for i in range(qr_flip.shape[0]):
		for k in range(qr_flip.shape[1]):
			if qr_flip[i][k] >= 0.001:
				x.append(i)
				z.append(k)

	scatter = plt.scatter(z, x, c='deepskyblue', s=5, marker='|')

	h1, _ = CS.legend_elements()
	h2, _ = CF.legend_elements()

	plt.legend([h1[0], scatter], [r"$q_c\quad[\frac{g}{kg}]$", r"$q_r\quad$"], loc='upper right')
	plt.savefig(f"graphs/qc+qr+th+u+w/qc+qr+th+u+w_{count}.png", dpi=300)
	plt.close()

	############################# qr+th+u+w ##########################
	# th = getData(np.loadtxt("outputs/th/th_" + str(t) + ".txt").tolist())
	# th_flip = th.swapaxes(0, 1)[:, 250:350+1]
	plt.figure(figsize=(10, 6))
	plt.tight_layout()
	plt.title(f"t = {t * dt} s", fontsize=14)
	plt.xlabel("x [km]")
	plt.ylabel("z [km]")
	plt.xticks(np.linspace(0, 100, 7), ["60", "65", "70", "75", "80", "85", "90"])
	plt.yticks([0, 10, 20, 30, 40, 50, 59], ["0", "2.5", "5", "7.5", "10", "12.5", "15"])
	CF = plt.contourf(th_flip, levels=np.arange(-16, 16+2,2), extend='both', cmap=cm.twilight_shifted)
	cbar = plt.colorbar(pad=0.05)
	cbar.set_ticks(np.arange(-16, 16+4, 4))
	cbar.set_label(r"$\theta^{'}$ [K]")

	X, Y = np.meshgrid(np.arange(101), np.arange(60))
	# u = getData(np.loadtxt("outputs/u/u_" + str(t) + ".txt").tolist())
	# u_flip = u.swapaxes(0, 1)[:, 250:350+1]
	# w = getData(np.loadtxt("outputs/w/w_" + str(t) + ".txt").tolist())
	# w_flip = w.swapaxes(0, 1)[:, 250:350+1]
	Q = plt.quiver(X[::2, ::2], Y[::2, ::2], u_flip[::2, ::2], w_flip[::2, ::2], angles='xy', units="width", scale=300)
	qk = plt.quiverkey(Q, 0.7, 0.9, 10, r'$10 \frac{m}{s}$', labelpos='E', coordinates='figure')

	# qr = getData2(np.loadtxt("outputs/qr/qr_" + str(t) + ".txt").tolist())
	# qr_flip = qr.swapaxes(0, 1)[:, 250:350+1]
	x, z = [], []
	for i in range(qr_flip.shape[0]):
		for k in range(qr_flip.shape[1]):
			if qr_flip[i][k] >= 0.001:
				x.append(i)
				z.append(k)

	scatter = plt.scatter(z, x, c='deepskyblue', s=5, marker='|')

	h2, _ = CF.legend_elements()

	plt.legend([scatter], [r"$q_r\quad$"], loc='upper right')
	plt.savefig(f"graphs/qr+th+u+w/qr+th+u+w_{count}.png", dpi=300)
	plt.close()

	############################# th ##########################
	# th = getData(np.loadtxt("outputs/th/th_" + str(t) + ".txt").tolist())
	# th_flip = th.swapaxes(0, 1)[:, 250:350+1]
	plt.figure(figsize=(10, 6))
	plt.tight_layout()
	plt.title(f"θ' [°C], t = {t * dt} s", fontsize=14)
	plt.xlabel("x [km]")
	plt.ylabel("z [km]")
	plt.xticks(np.linspace(0, 100, 7), ["60", "65", "70", "75", "80", "85", "90"])
	plt.yticks([0, 10, 20, 30, 40, 50, 59], ["0", "2.5", "5", "7.5", "10", "12.5", "15"])
	plt.contourf(th_flip, levels=np.arange(-16, 16+2,2), extend='both', cmap=cm.twilight_shifted)
	cbar = plt.colorbar(pad=0.05)
	cbar.set_ticks(np.arange(-16, 16+2,2))
	plt.savefig(f"graphs/th/th_{count}.png", dpi=300)
	plt.close()

	############################# zeta ##########################
	zeta = np.loadtxt("outputs/zeta/zeta_" + str(t) + ".txt").reshape([nx, nz], order='F')
	zeta_flip = zeta.swapaxes(0, 1)[:, 250:350+1]
	plt.figure(figsize=(10, 6))
	plt.tight_layout()
	plt.title(r"$\zeta$ $[\frac{1}{s}]$" f",  t = {t * dt} s", fontsize=14)
	plt.xlabel("x [km]")
	plt.ylabel("z [km]")
	plt.xticks([0, 100, 200, 300, 400, 500, 599], ["0", "25", "50", "75", "100", "125", "150"])
	plt.yticks([0, 10, 20, 30, 40, 50, 59], ["0", "2.5", "5", "7.5", "10", "12.5", "15"])
	plt.contourf(zeta_flip, levels=np.arange(-0.2, 0.2+0.02, 0.02), extend='both', cmap=cm.twilight_shifted)
	cbar = plt.colorbar(pad=0.05)
	cbar.set_ticks(np.arange(-0.2, 0.2+0.02, 0.02))
	plt.savefig(f"graphs/zeta/zeta_{count}.png", dpi=300)
	plt.close()



	############################# u ##########################
	# u = getData(np.loadtxt("outputs/u/u_" + str(t) + ".txt").tolist())
	# u_flip = u.swapaxes(0, 1)[:, 250:350+1]
	plt.figure(figsize=(10, 6))
	plt.tight_layout()
	plt.title(r"u $[\frac{m}{s}]$" f",  t = {t * dt} s", fontsize=14)
	plt.xlabel("x [km]")
	plt.ylabel("z [km]")
	plt.xticks(np.linspace(0, 100, 7), ["60", "65", "70", "75", "80", "85", "90"])
	plt.yticks([0, 10, 20, 30, 40, 50, 59], ["0", "2.5", "5", "7.5", "10", "12.5", "15"])
	plt.contourf(u_flip, levels=np.arange(-36, 36+4, 4), extend='both', cmap=cm.twilight_shifted)
	cbar = plt.colorbar(pad=0.05)
	cbar.set_ticks(np.arange(-36, 36+4, 4))
	plt.savefig(f"graphs/u/u_{count}.png", dpi=300)
	plt.close()


	############################# w ##########################
	# w = getData(np.loadtxt("outputs/w/w_" + str(t) + ".txt").tolist())
	# w_flip = w.swapaxes(0, 1)[:, 250:350+1]
	plt.figure(figsize=(10, 6))
	plt.tight_layout()
	plt.title(r"w $[\frac{m}{s}]$" f",  t = {t * dt} s", fontsize=14)
	plt.xlabel("x [km]")
	plt.ylabel("z [km]")
	plt.xticks([0, 100, 200, 300, 400, 500, 599], ["0", "25", "50", "75", "100", "125", "150"])
	plt.yticks([0, 10, 20, 30, 40, 50, 59], ["0", "2.5", "5", "7.5", "10", "12.5", "15"])
	plt.contourf(w_flip, levels=np.arange(-36, 36+4, 4), extend='both', cmap=cm.twilight_shifted)
	cbar = plt.colorbar(pad=0.05)
	cbar.set_ticks(np.arange(-36, 36+4, 4))
	plt.savefig(f"graphs/w/w_{count}.png", dpi=300)
	plt.close()


	############################# qv ##########################
	qv = np.loadtxt("outputs/qv/qv_" + str(t) + ".txt").reshape([nx, nz], order='F')
	qv_flip = qv.swapaxes(0, 1)[:, 250:350+1]
	plt.figure(figsize=(10, 6))
	plt.tight_layout()
	plt.title(r"$q_v$ $[\frac{kg}{kg}]$" + f",  t = {t * dt} s", fontsize=14)
	plt.xlabel("x [km]")
	plt.ylabel("z [km]")
	plt.xticks(np.linspace(0, 100, 7), ["60", "65", "70", "75", "80", "85", "90"])
	plt.yticks([0, 10, 20, 30, 40, 50, 59], ["0", "2.5", "5", "7.5", "10", "12.5", "15"])
	plt.contourf(qv_flip, extend='max', levels=np.linspace(0, 0.014, 7))
	cbar = plt.colorbar(pad=0.05)
	cbar.set_ticks(np.linspace(0, 0.014, 7))
	plt.savefig(f"graphs/qv/qv_{count}.png", dpi=300)
	plt.close()

	############################# qc ##########################
	# qc = getData2(np.loadtxt("outputs/qc/qc_" + str(t) + ".txt").tolist())
	# qc_flip = qc.swapaxes(0, 1)[:, 250:350+1]
	plt.figure(figsize=(10, 6))
	plt.tight_layout()
	plt.title(r"$q_c$ $[\frac{kg}{kg}]$" + f",  t = {t * dt} s", fontsize=14)
	plt.xlabel("x [km]")
	plt.ylabel("z [km]")
	plt.xticks(np.linspace(0, 100, 7), ["60", "65", "70", "75", "80", "85", "90"])
	plt.yticks([0, 10, 20, 30, 40, 50, 59], ["0", "2.5", "5", "7.5", "10", "12.5", "15"])
	plt.contourf(qc_flip, extend='max', levels=np.linspace(0, 0.006, 7))
	cbar = plt.colorbar(pad=0.05)
	cbar.set_ticks(np.linspace(0, 0.006, 7))
	plt.savefig(f"graphs/qc/qc_{count}.png", dpi=300)
	plt.close()

	############################# qr ##########################
	# qr = getData2(np.loadtxt("outputs/qr/qr_" + str(t) + ".txt").tolist())
	# qr_flip = qr.swapaxes(0, 1)[:, 250:350+1]
	plt.figure(figsize=(10, 6))
	plt.tight_layout()
	plt.title(r"$q_r$ $[\frac{kg}{kg}]$" + f",  t = {t * dt} s", fontsize=14)
	plt.xlabel("x [km]")
	plt.ylabel("z [km]")
	plt.xticks(np.linspace(0, 100, 7), ["60", "65", "70", "75", "80", "85", "90"])
	plt.yticks([0, 10, 20, 30, 40, 50, 59], ["0", "2.5", "5", "7.5", "10", "12.5", "15"])
	plt.contourf(qr_flip, extend='max', levels=np.linspace(0, 0.04, 9))
	cbar = plt.colorbar(pad=0.05)
	cbar.set_ticks(np.linspace(0, 0.04, 9))
	plt.savefig(f"graphs/qr/qr_{count}.png", dpi=300)
	plt.close()

	############################# qc+qr ##########################
	# qr = np.loadtxt("outputs/qr/qr_" + str(t) + ".txt").reshape([nx, nz], order='F')
	qr_flip = qr.swapaxes(0, 1)[:, 250:350+1]
	plt.figure(figsize=(10, 6))
	plt.tight_layout()
	plt.title(r"$q_r$ $[\frac{kg}{kg}]$ (shading)" + r" & $q_{c}$ $[\frac{kg}{kg}]$ (contour)" + f",  t = {t * dt} s", fontsize=14)
	plt.xlabel("x [km]")
	plt.ylabel("z [km]")
	plt.xticks(np.linspace(0, 100, 7), ["60", "65", "70", "75", "80", "85", "90"])
	plt.yticks([0, 10, 20, 30, 40, 50, 59], ["0", "2.5", "5", "7.5", "10", "12.5", "15"])
	plt.contourf(qr_flip, extend='max', levels=np.linspace(0, 0.04, 9))
	cbar = plt.colorbar(pad=0.05)
	cbar.set_ticks(np.linspace(0, 0.04, 9))

	# qc = np.loadtxt("outputs/qc/qc_" + str(t) + ".txt").reshape([nx, nz], order='F')
	qc_flip = qc.swapaxes(0, 1)[:, 250:350+1]
	CS = plt.contour(qc_flip, extend='max', levels=[0.0000001, 0.001, 0.002, 0.003, 0.004, 0.005], colors='red', linewidths=1.25)
	plt.clabel(CS, inline=1, fontsize=8, fmt="%.4f")
	plt.savefig(f"graphs/qc+qr/qc+qr_{count}.png", dpi=300)
	plt.close()

	############################# qv+qc ##########################
	# qv = getData(np.loadtxt("outputs/qv/qv_" + str(t) + ".txt").tolist())
	# qv_flip = qv.swapaxes(0, 1)[:, 250:350+1]
	plt.figure(figsize=(10, 6))
	plt.tight_layout()
	plt.title(r"$q_v$ $[\frac{kg}{kg}]$" + r" & $q_{c}$ $[\frac{kg}{kg}]$" + f",  t = {t * dt} s", fontsize=14)
	plt.xlabel("x [km]")
	plt.ylabel("z [km]")
	plt.xticks(np.linspace(0, 100, 7), ["60", "65", "70", "75", "80", "85", "90"])
	plt.yticks([0, 10, 20, 30, 40, 50, 59], ["0", "2.5", "5", "7.5", "10", "12.5", "15"])
	plt.contourf(qv_flip, extend='max', levels=np.linspace(0, 0.014, 7))
	cbar = plt.colorbar(pad=0.05)
	cbar.set_ticks(np.linspace(0, 0.014, 7))

	# qc = getData2(np.loadtxt("outputs/qc/qc_" + str(t) + ".txt").tolist())
	# qc_flip = qc.swapaxes(0, 1)[:, 250:350+1]
	CS = plt.contour(qc_flip, extend='max', levels=np.linspace(0, 0.003, 4), colors='red', linewidths=1.25)
	plt.clabel(CS, inline=1, fontsize=8, fmt="%.4f")
	plt.savefig(f"graphs/qv+qc/qv+qc_{count}.png", dpi=300)
	plt.close()
	

	count += 1
	print(t)



"""
for t in range(0, 20001, 250):
	th = getData(np.loadtxt("outputs/th/th_" + str(t) + ".txt").tolist())
	th_flip = th.swapaxes(0, 1)[:, 50:]
	# TH = PlotMethod1D(th_flip)
	plt.tight_layout()
	plt.figure(figsize=(10, 6))
	plt.title(f"θ' [°C], t = {t * dt}", fontsize=14)
	plt.xlabel("x [km]")
	plt.ylabel("z [km]")
	plt.xticks([0, 100, 200, 300, 400, 500, 599], ["0", "25", "50", "75", "100", "125", "150"])
	plt.yticks([0, 10, 20, 30, 40, 50, 59], ["0", "2.5", "5", "7.5", "10", "12.5", "15"])
	plt.contourf(th_flip, levels=np.arange(-32, 32+2,2), extend='both')
	cbar = plt.colorbar(pad=0.05)
	cbar.set_ticks(np.arange(-32, 32+4, 4))
	plt.savefig(f"graphs/th/th_{count}.png")
	plt.close()
	# TH.contourf(title=f"θ', t = {t * dt}", show=False, savename=f"graphs/th/th_{count}.png", levels=np.arange(-30, 30+2,2))

	# ADV U
	# th2d = np.asarray(th_flip[10, :])
	# TH2d = PlotMethod2D(np.arange(128), th2d)
	# TH2d.plot2D(show=False, savename=f"graphs/th2d/th_{count}.png")

	# ADV W
	# th2d = np.asarray(th_flip[:, int(nx / 2) - 1])
	# TH2d = PlotMethod2D(np.arange(60), th2d)
	# TH2d.plot2D(show=False, savename=f"graphs/th2d/th_{count}.png")


	zeta = getData(np.loadtxt("outputs/zeta/zeta_" + str(t) + ".txt").tolist())
	zeta_flip = zeta.swapaxes(0, 1)
	ZETA = PlotMethod1D(zeta_flip)
	plt.tight_layout()
	plt.figure(figsize=(10, 6))
	plt.title(r"$\zeta$ $[\frac{1}{s}]$" f",  t = {t * dt}", fontsize=14)
	plt.xlabel("x [km]")
	plt.ylabel("z [km]")
	plt.xticks([0, 100, 200, 300, 400, 500, 599], ["0", "25", "50", "75", "100", "125", "150"])
	plt.yticks([0, 10, 20, 30, 40, 50, 59], ["0", "2.5", "5", "7.5", "10", "12.5", "15"])
	plt.contourf(zeta_flip, levels=np.arange(-0.2, 0.2+0.02, 0.02), extend='both')
	cbar = plt.colorbar(pad=0.05)
	cbar.set_ticks(np.arange(-0.2, 0.2+0.02, 0.02))
	plt.savefig(f"graphs/zeta/zeta_{count}.png")
	plt.close()

	u = getData(np.loadtxt("outputs/u/u_" + str(t) + ".txt").tolist())
	u_flip = u.swapaxes(0, 1)
	U = PlotMethod1D(u_flip)
	# U.contourf(figsize=(10, 6), title=f"u, t = {t * dt}", show=False, savename=f"graphs/u/u_{count}.png")
	plt.tight_layout()
	plt.figure(figsize=(10, 6))
	plt.title(r"u $[\frac{m}{s}]$" f",  t = {t * dt}", fontsize=14)
	plt.xlabel("x [km]")
	plt.ylabel("z [km]")
	plt.xticks([0, 100, 200, 300, 400, 500, 599], ["0", "25", "50", "75", "100", "125", "150"])
	plt.yticks([0, 10, 20, 30, 40, 50, 59], ["0", "2.5", "5", "7.5", "10", "12.5", "15"])
	plt.contourf(u_flip, levels=np.arange(-36, 36+4, 4), extend='both')
	cbar = plt.colorbar(pad=0.05)
	cbar.set_ticks(np.arange(-36, 36+4, 4))
	plt.savefig(f"graphs/u/u_{count}.png")
	plt.close()

	w = getData(np.loadtxt("outputs/w/w_" + str(t) + ".txt").tolist())
	w_flip = w.swapaxes(0, 1)
	W = PlotMethod1D(w_flip)
	# W.contourf(figsize=(10, 6), title=f"w, t = {t * dt}", show=False, savename=f"graphs/w/w_{count}.png")
	plt.tight_layout()
	plt.figure(figsize=(10, 6))
	plt.title(r"w $[\frac{m}{s}]$" f",  t = {t * dt}", fontsize=14)
	plt.xlabel("x [km]")
	plt.ylabel("z [km]")
	plt.xticks([0, 100, 200, 300, 400, 500, 599], ["0", "25", "50", "75", "100", "125", "150"])
	plt.yticks([0, 10, 20, 30, 40, 50, 59], ["0", "2.5", "5", "7.5", "10", "12.5", "15"])
	plt.contourf(w_flip, levels=np.arange(-36, 36+4, 4), extend='both')
	cbar = plt.colorbar(pad=0.05)
	cbar.set_ticks(np.arange(-36, 36+4, 4))
	plt.savefig(f"graphs/w/w_{count}.png")
	plt.close()

	qv = getData(np.loadtxt("outputs/qv/qv_" + str(t) + ".txt").tolist())
	qv_flip = qv.swapaxes(0, 1)
	QV = PlotMethod1D(qv_flip)
	# QV.contourf(figsize=(10, 6), title=r"$q_v$ $[\frac{kg}{kg}]$" + f",  t = {t * dt}", show=False, 
	# 			savename=f"graphs/qv/qv_{count}.png", levels=np.arange(-0.006, 0.006+0.0002, 0.0002))
	plt.tight_layout()
	plt.figure(figsize=(10, 6))
	plt.title(r"$q_v$ $[\frac{kg}{kg}]$" + f",  t = {t * dt}", fontsize=14)
	plt.xlabel("x [km]")
	plt.ylabel("z [km]")
	plt.xticks([0, 100, 200, 300, 400, 500, 599], ["0", "25", "50", "75", "100", "125", "150"])
	plt.yticks([0, 10, 20, 30, 40, 50, 59], ["0", "2.5", "5", "7.5", "10", "12.5", "15"])
	plt.contourf(qv_flip, extend='both', levels=np.arange(-0.006, 0.006+0.0002, 0.0002))
	cbar = plt.colorbar(pad=0.05)
	cbar.set_ticks(np.arange(-0.006, 0.006+0.0008, 0.0008))
	plt.savefig(f"graphs/qv/qv_{count}.png")
	plt.close()

	qc = getData2(np.loadtxt("outputs/qc/qc_" + str(t) + ".txt").tolist())
	qc_flip = qc.swapaxes(0, 1)
	QC = PlotMethod1D(qc_flip)
	# QC.contourf(figsize=(10, 6), title=f"qc, t = {t * dt}", show=False, savename=f"graphs/qc/qc_{count}.png")
	plt.tight_layout()
	plt.figure(figsize=(10, 6))
	plt.title(r"$q_c$ $[\frac{kg}{kg}]$" + f",  t = {t * dt}", fontsize=14)
	plt.xlabel("x [km]")
	plt.ylabel("z [km]")
	plt.xticks([0, 100, 200, 300, 400, 500, 599], ["0", "25", "50", "75", "100", "125", "150"])
	plt.yticks([0, 10, 20, 30, 40, 50, 59], ["0", "2.5", "5", "7.5", "10", "12.5", "15"])
	plt.contourf(qc_flip, extend='both', levels=np.arange(0, 0.006+0.0002, 0.0002))
	cbar = plt.colorbar(pad=0.05)
	cbar.set_ticks(np.arange(0, 0.006+0.0004, 0.0004))
	plt.savefig(f"graphs/qc/qc_{count}.png")
	plt.close()

	qr = getData2(np.loadtxt("outputs/qr/qr_" + str(t) + ".txt").tolist())
	qr_flip = qr.swapaxes(0, 1)
	QR = PlotMethod1D(qr_flip)
	# QR.contourf(figsize=(10, 6), title=f"qr, t = {t * dt}", show=False, savename=f"graphs/qr/qr_{count}.png", levels=np.arange(0, 0.08+0.002, 0.002))
	plt.tight_layout()
	plt.figure(figsize=(10, 6))
	plt.title(r"$q_r$ $[\frac{kg}{kg}]$" + f",  t = {t * dt}", fontsize=14)
	plt.xlabel("x [km]")
	plt.ylabel("z [km]")
	plt.xticks([0, 100, 200, 300, 400, 500, 599], ["0", "25", "50", "75", "100", "125", "150"])
	plt.yticks([0, 10, 20, 30, 40, 50, 59], ["0", "2.5", "5", "7.5", "10", "12.5", "15"])
	plt.contourf(qr_flip, extend='both', levels=np.arange(0, 0.08+0.002, 0.002))
	cbar = plt.colorbar(pad=0.05)
	cbar.set_ticks(np.arange(0, 0.08+0.004, 0.004))
	plt.savefig(f"graphs/qr/qr_{count}.png")
	plt.close()


	count += 1
	print(t)

"""

