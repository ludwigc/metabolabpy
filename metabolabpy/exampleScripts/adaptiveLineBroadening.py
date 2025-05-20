message = True
s = self.nd.s
min_lb = 0.3
for k in range(len(self.nd.nmrdat[s])):
    self.nd.nmrdat[s][k].proc.window_type[0] = 1
    self.nd.nmrdat[s][k].proc.lb[0] = min_lb

self.ft_all()
if message == True:
    print('\n--------------------------------------------------\n')
msg = self.fit_tmsp_all(message)
if message == True:
    print('\n--------------------------------------------------\n')
lw = []
for k in range(len(self.nd.nmrdat[s])):
    lw.append(self.nd.nmrdat[s][k].tmsp_linewidth)

lw = np.array(lw)
idx = np.where(lw == np.max(lw))[0]
max_lw = lw[idx]
for k in range(len(self.nd.nmrdat[s])):
    self.nd.nmrdat[s][k].proc.lb[0] = min_lb + max_lw - lw[k]

self.ft_all()
self.fit_tmsp_all()
if message == True:
    print('\n--------------------------------------------------\n')

self.plot_spc()
self.show_console()