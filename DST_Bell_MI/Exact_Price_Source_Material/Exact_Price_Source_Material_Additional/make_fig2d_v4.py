"""Regenerates panel (d) of Fig. 2 (arXiv:2608.18886 v4) from the ten exact-profile
harmonics a_k (as in protocol_sensitivity_v2_1.py), normalized to the mean, and pastes it
over the lower-right quadrant of the v3 figure. Usage: python make_fig2d_v4.py paper_fig1_mechanics.png
Writes panel_d.png and paper_fig1_mechanics_v4.png. v3's panel (d) was normalized to the
CHSH-angle value, which made the fringe read as zero there."""
import numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'mathtext.fontset':'dejavusans'})
from PIL import Image
A_K=np.array([5.773052e-02,9.834086e-03,3.100640e-03,1.331055e-03,6.853649e-04,
              3.973659e-04,2.503605e-04,1.677163e-04,1.177647e-04,8.582444e-05])
k=np.arange(1,11)
R=lambda d: 1+np.sum(A_K[:,None]*np.cos(4*k[:,None]*np.deg2rad(d)[None,:]),axis=0)
D=np.linspace(0,90,901)
G='#1d9e75'; O='#d85a30'; C='#555577'
W,H,dpi=977,702,170
fig=plt.figure(figsize=(W/dpi,H/dpi),dpi=dpi)
ax=fig.add_axes([122.5/W,1-604.5/H,(950.5-122.5)/W,(599.5-60.5)/H])
ax.plot(D,R(D),color=G,lw=2,label='this class: cos 4$\\Delta$ fringe')
ax.axhline(1.0,color=O,ls='--',lw=1.5,label='quantum mechanics: flat')
rc=R(np.array([22.5]))[0]
for x in (22.5,67.5):
    ax.axvline(x,color=C,ls=':',lw=1)
    ax.plot([x],[rc],'o',color=C,ms=4,zorder=5)
ax.plot([22.5,67.5],[rc,rc],color=C,ls=':',lw=1)
ax.text(24,1.055,'CHSH angle:\n%.3f of mean'%rc,color=C,fontsize=8,va='center')
ax.text(52.5,1.055,'CHSH angle:\n%.3f of mean'%rc,color=C,fontsize=8,va='center')
ax.text(45,rc-0.004,'equal, not zero',color=C,fontsize=8,ha='center',va='top')
r0,r45=R(np.array([0.0]))[0],R(np.array([45.0]))[0]
ax.plot([0,45],[r0,r45],'s',mfc='white',mec=G,ms=5,zorder=6)
ax.text(45,r45-0.003,'0° vs 45°: ratio %.2f'%(r0/r45),color=G,fontsize=8,ha='center',va='top')
ax.set_xlim(-4.5,94.5); ax.set_ylim(0.935,1.09)
ax.set_xlabel('relative angle $\\Delta$ (deg)')
ax.set_ylabel('total coincidence rate (norm. to mean)')
ax.set_title('(d) The fourfold-sum fringe: equal values,\nnot zero, at the two CHSH relative angles',fontsize=10)
ax.legend(loc='upper center',fontsize=8)
fig.savefig('panel_d.png',dpi=dpi)
print('R0',r0,'R45',r45,'R22.5',rc,'ratio',r0/r45)
import sys
src=sys.argv[1] if len(sys.argv)>1 else 'paper_fig1_mechanics.png'
orig=Image.open(src).convert('RGBA')
p=Image.open('panel_d.png').convert('RGBA'); print(p.size)
new=orig.copy(); new.paste(p,(978,760))
new.save('paper_fig1_mechanics_v4.png',dpi=(170,170))
