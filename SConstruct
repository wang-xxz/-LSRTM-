# —*- coding=utf-8 -*-
from rsf.proj import *
from rsf.recipes import fdmod


def igrey3d(par):
    return '''
    byte allpos=n gainpanel=all |
    grey3 flat=n frame1=2000 frame2=200 frame3=10
    title= titlesz=10 font=3
    label2=Trace unit2=km lable3=Shot unit3= point1=0.85 point2=0.85
    labelsz=5
    %s
    ''' % par


def igrey2d(par):
    'grey 2d image for RTM images'
    return '''
    grey title=
    scalebar=y bartype=v barwidth=0.2
    screenratio=0.7 screenht=6
    labelsz=5 labelfat=2
    %s
    ''' % par


def mgrey2d(par):
    'grey 2d image for encoded and cross-talk matrix images'
    return '''
    grey gainpanel=a bias=0.5
    screenht=10 font=3
    wantscalebar=n parallel2=n labelsz=5 labelfat=3 titlesz=10 titlefat=3
    %s
    ''' % par


par_ = dict(nx=2301, ox=0, dx=0.004, lx='x', ux='km',
            nz=751, oz=0, dz=0.004, lz='z', uz='km',
            )

# Download velocity model
#Fetch('marmvel.hh', 'marm', usedatapath=False)
Flow('marm-vel', 'marmvel.hh',
     '''
    dd form=native | put d1=%(dz)f d2=%(dx)f
    label1=%(lz)s  unit1=%(uz)s label2=%(lx)s unit2=%(ux)s title= |
    scale rscale=0.001 |
    window j1=2 j2=2
    out=stdout
    ''' % par_)

# 将模型拓展20个网格点
# Flow('marm_vel', 'zvel',
#      '''
#     window n1=20 | math output=1.5 | cat axis=1 ${SOURCE}
#     ''')

# parameters dict
par = dict(
    nx=1151, ox=0, dx=0.008, lx='x', ux='km',
    nz=376,  oz=0, dz=0.008, lz='z', uz='km',
    nt=5000, ot=0, dt=0.0008, lt='t', ut='s',
    kt=100,
    nb=30,
    freq=15
)

# extend model
rtmpar = dict(
    lft=0,
    rht=0,
    top=20,
    bot=20
)
par['nx'] = par['nx'] + rtmpar['lft'] + rtmpar['rht']
par['nz'] = par['nz'] + rtmpar['top'] + rtmpar['bot']
par['ox'] = par['ox'] - rtmpar['lft'] * par['dx']
par['oz'] = par['oz'] - rtmpar['top'] * par['dz']

fdmod.param(par)

winpar = 'min1=0 max1=3 min2=0 n1=376 n2=1151'

Flow('marm-vel-ext', 'marm-vel',
     'pad2 left=%(lft)d right=%(rht)d top=%(top)d bottom=%(bot)d' % rtmpar)
Plot('marm-vel-ext', fdmod.cgrey('color=j allpos=y bias=1.5 title=', par))

Result('marm-vel', 'marm-vel-ext',
       'window %s |' % winpar +
       igrey2d('color=i mean=y wantaxis2=y scalebar=y bartype=v \
        barwidth=0.2 barunit=km/s barlabel=Velocity'))

# creat constant velocity for substract direct wave
cut = 2
Flow('marm-vel-con', 'marm-vel-ext', 'window n1=%d | pad2 bottom=%d'
     % (cut, par['nz'] - cut))


# Smoothed velocity
Flow('marm-vel-ext-smo', 'marm-vel-ext', 'smooth rect1=10 rect2=10 repeat=2')

# Density
Flow('marm-dens-ext', 'marm-vel-ext', 'math output=1')

# compute reflectivity
Flow('marm-ref', 'marm-vel',
     '''
    depth2time velocity=$SOURCE nt=%(nt)d dt=%(dt)f |
    ai2refl |
    ricker1 frequency=%(freq)d |
    time2depth velocity=$SOURCE
    ''' % par)
Result('marm-ref', igrey2d(
    'title= wantscalebar=y minval=-0.02 maxval=0.02 pclip=99\
    barlabel=Amplitude') + winpar)


#
# source Wavelet
# Flow('marm_wlt', None,
#      '''
#      spike n1=%(nt)d o1=0 d1=%(dt)f k1=400 mag=1 nsp=1 |
#      ricker1 frequency=15 |
#      window n1=%(nt)d |
#      scale axis=123 |
#      put label1=%(lt)s unit1=%(ut)s
#      ''' % par)
#
fdmod.wavelet('wlt_', par['freq'], par)
Flow('marm-wlt', 'wlt_', 'transp')

Result('marm-wlt', 'window n2=2000 |' + fdmod.waveplot('', par))
# =============================================================================
# Conventional RTM result(FD and born data)
# for comparetion
#
# 3. 炮点和接收点坐标设置
# ======================================================================
# 水平接收点坐标(z=0)
jrx = 1
par['drx'] = jrx * par['dx']
par['orx'] = par['ox']    # origin
par['nrx'] = par['nx'] / jrx   # number
par['rz'] = par['oz'] + 0.025

Flow('rx', None, 'math n1=%(nrx)d d1=%(drx)f o1=%(orx)f output=x1' % par)
Flow('rz', None, 'math n1=%(nrx)d d1=%(drx)f o1=%(orx)f output=%(rz)f' % par)
Flow('rxz', 'rx rz', 'cat axis=2 space=n ${SOURCES[1]} | transp ')
Plot('rxz', fdmod.rrplot('', par))
Result('rxz', ['marm-vel-ext', 'rxz'], 'Overlay')

# 水平炮点坐标
par['jsx'] = 18
par['dsx'] = par['dx'] * par['jsx']
par['osx'] = par['ox'] + 0.1
par['nsx'] = 60
par['sz'] = par['oz'] + 0.01

Flow('sx', None, 'math n1=%(nsx)d d1=%(dsx)f o1=%(osx)f output=x1' % par)
Flow('sz', None, 'math n1=%(nsx)d d1=%(dsx)f o1=%(osx)f output=%(sz)f' % par)
Flow('sxz', 'sx sz', 'cat axis=2 space=n ${SOURCES[1]} | transp')
Plot('sxz', fdmod.ssplot('', par))
Result('sxz', ['marm-vel-ext', 'sxz'], 'Overlay')

# 4. RTM单炮成像
# ======================================================================
# prog = Program(Split('Mprertm2d.c prertm2d.c'), PROGSUFFIX='.x')
# exe = str(prog[0])

# conventional RTM program
prog1 = Program(Split('../SRC/Mprertm2d_v02.c ../SRC/prertm2d_v02.c'),
                PROGSUFFIX='.x')
exe1 = str(prog1[0])
# LSRTM with regularization
prog3 = Program(Split('../SRC/Mlsrtmsr.c ../SRC/lsrtmsr.c \
                ../SRC/laplac2.c'),
                PROGSUFFIX='.x',
                LIBS=['rsfpwd', 'rsf', 'm', 'gomp', 'cblas'])
exe3 = str(prog3[0])
# ==================================================

# Flow('marm-img', ['marm-data', 'marm-wlt', 'marm-vel-ext-smo',
#                   'sxz', 'rxz', exe1],
#      '''
#     ${SOURCES[5].abspath} verb=y adj=y nb=100
#     wlt=${SOURCES[1]} vel=${SOURCES[2]} sou=${SOURCES[3]} rec=${SOURCES[4]}
#     ''')
# Result('marm-img',
#        'window %s | laplac | igrad |' % winpar +
#        igrey2d('color=i scalebar=n pclip=99'))

Flow('marm-img-illum', ['marm-data', 'marm-wlt', 'marm-vel-ext-smo',
                        'sxz', 'rxz', exe1],
     '''
    ${SOURCES[5].abspath} verb=y adj=y nb=100
    wlt=${SOURCES[1]} vel=${SOURCES[2]} sou=${SOURCES[3]} rec=${SOURCES[4]}
    ''')
Result('marm-img-illum',
       'window %s | laplac | math output=1e4*input|' % winpar +
       igrey2d('color=i wantscalebar=y pclip=99 minval=-0.05 maxval=0.05\
       barlabel=Amplitude barunit="*10e-4"'))

Result('marm-img-wgt', 'marm-img-illum',
       'window %s | laplac | smooth rect1=5 rect2=1 | bandpass flo=5|' % winpar +
       igrey2d('color=i wantscalebar=y pclip=99'))

Result('marm-img-win', 'marm-img-illum',
       'window min1=0.2 max1=1.5 min2=4 max2=7 | laplac | smooth rect1=5 rect2=1 | bandpass flo=5|' +
       igrey2d('color=i scalebar=n pclip=99'))
# ------------------
Flow(['marm-invs-img', 'marm-invs-img-err'],
     ['marm-data',
      'marm-wlt',
      'marm-vel-ext-smo',
      'sxz',
      'rxz',
      exe3],
     '''
     ${SOURCES[5].abspath} verb=y nb=100 nss=%d niter=5
     wlt=${SOURCES[1]} vel=${SOURCES[2]} sou=${SOURCES[3]} rec=${SOURCES[4]}
     error=${TARGETS[1]}
     ''' % par['nsx'])
Result('marm-invs-img',
       'window %s | ' % winpar +
       igrey2d('color=i wantscalebar=y pclip=99'))


End()
