"""Regenerate the VLEACH 2.2a reference cases (data/reference/vleach/2.2a/) with EPA's executable.

Not part of the test run. Requirements: Wine (wine32), msitools (msiextract).
Steps:
  1. extract VleachForGUI.exe (Intel Visual Fortran console program) and the EPA sample
     problem from docs/references/vleach_2.2a_setup.exe: the setup is an InstallShield
     package; its overlay holds 4 files (0x0409.ini, MSIEng.isc, 'VLEACH Model Interface.msi',
     Setup.ini), each XOR-encoded with a key derived from its file name (function
     extract_setup below); the .msi is unpacked with msiextract;
  2. write each case as XXXXXX.INP (the fixed input name of the executable; format:
     VLEACH UG section 6, Figure 4), run the executable, and keep the input and the
     .OUT, .PRF, .PRM and GWIMP.DAT files (line ends converted to LF).
The executable reproduces the published Sample.out byte for byte (checked).

Run:  python modules/vadose/tests/generate_vleach_reference.py --workdir /tmp/vleach
"""
import argparse, os, shutil, struct, subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
SETUP = ROOT / "docs/references/vleach_2.2a_setup.exe"
OUT = ROOT / "data/reference/vleach/2.2a"


def extract_setup(setup: Path, dest: Path) -> Path:
    """Unpack the InstallShield overlay and the MSI; return the folder 'VLEACH Model'."""
    import pefile
    d = setup.read_bytes()
    off = pefile.PE(str(setup)).get_overlay_data_start_offset()
    assert d[off:off + 14] == b"InstallShield\x00"
    nfiles = struct.unpack_from("<H", d, off + 14)[0]
    p, magic = off + 0x2E, (0x13, 0x35, 0x86, 0x07)
    dest.mkdir(parents=True, exist_ok=True)
    for _ in range(nfiles):
        name = d[p:p + 260].split(b"\0")[0].decode("latin-1")
        flags, _x, flen = struct.unpack_from("<III", d, p + 260)
        p += 312
        data = d[p:p + flen]; p += flen
        key = [ord(ch) ^ magic[i % 4] for i, ch in enumerate(name)]
        if flags & 2:
            data = bytes((~(key[i % len(key)] ^ (((x << 4) | (x >> 4)) & 0xFF))) & 0xFF for i, x in enumerate(data))
        (dest / name).write_bytes(data)
    subprocess.run(["msiextract", "-C", str(dest), str(dest / "VLEACH Model Interface.msi")], check=True,
                   capture_output=True)
    return dest / "Program Files" / "VLEACH Model"


def g(x): return f"{x:>10}"
def inp(title, delt, stime, ptime, prtime, koc, kh, cmax, dair, polys):
    L=[title, f"{len(polys):3d}", g(delt)+g(stime)+g(ptime)+g(prtime), g(koc)+g(kh)+g(cmax)+g(dair)]
    for p in polys:
        L.append(p['title'])
        L.append(''.join(g(p[k]) for k in ('area','delz','q','rhob','por','theta','foc')))
        L.append(g(p['cinf'])+g(p['catm'])+g(p['cgw']))
        L.append(f"{p['ncell']:5d}Y{g(stime)}")
        for j1,j2,c in p['layers']: L.append(f"{j1:5d}{j2:5d}{g(c)}")
    return '\r\n'.join(L)+'\r\n'
def P(title,ncell,layers,**kw):
    d=dict(title=title,area=1000,delz=1,q=1,rhob=1.6,por=0.4,theta=0.3,foc=0.005,cinf=0,catm=0,cgw=0,ncell=ncell,layers=layers); d.update(kw); return d
CASES={
 'case1_boundaries': inp('Case 1: TCE, recharge conc., atmosphere and groundwater gas boundaries',1,50,5,10,100,0.4,1100,0.7,
     [P('Column',30,[(1,10,500),(11,20,100),(21,30,0)],cinf=0.05,catm=0.001,cgw=0.01)]),
 'case2_impermeable': inp('Case 2: PCE, both boundaries impermeable, fine cells',0.1,20,2,5,265,1.0,200,0.62,
     [P('Fine',40,[(1,5,2000),(6,40,0)],delz=0.5,q=0.5,catm=-1,cgw=-1,foc=0.002,rhob=1.7,por=0.35,theta=0.15)]),
 'case3_three_polygons': inp('Case 3: benzene, three polygons',0.5,30,5,15,60,0.22,1790,0.76,
     [P('Sand',20,[(1,8,800),(9,20,0)],area=5000,q=2,por=0.38,theta=0.12,foc=0.001,rhob=1.65),
      P('Loam',20,[(1,4,1500),(5,20,0)],area=2500,q=0.8,por=0.45,theta=0.3,foc=0.01,rhob=1.45),
      P('Clay',10,[(1,10,300)],area=800,q=0.2,por=0.5,theta=0.42,foc=0.02,rhob=1.3,delz=2)]),
 'case4_diffusion_only': inp('Case 4: vinyl chloride, no recharge (diffusion only)',0.05,5,0.5,1,22,1.1,8800,0.92,
     [P('Dry',25,[(1,9,0),(10,15,1000),(16,25,0)],q=0,theta=0.08,por=0.36,foc=0.001)]),
 'case5_advection_only': inp('Case 5: non-volatile (H = 0, Dair = 0), advection only',0.2,40,4,8,30000,0.0,10,0.0,
     [P('Column',50,[(1,5,5000),(6,50,0)],q=3,foc=0.001,catm=-1,cgw=-1)]),
}
ap = argparse.ArgumentParser()
ap.add_argument("--workdir", default="/tmp/vleach")
work = Path(ap.parse_args().workdir)
model = extract_setup(SETUP, work / "setup")
EXE = model / "SystemFiles" / "VleachForGUI.exe"
env = dict(os.environ, WINEDEBUG="-all", WINEPREFIX=str(work / "wineprefix"))
for name,text in CASES.items():
    d=f'{OUT}/{name}'; shutil.rmtree(d,ignore_errors=True); os.makedirs(d)
    wd=str(work / 'run'); shutil.rmtree(wd,ignore_errors=True); os.makedirs(wd); shutil.copy(EXE,wd)
    open(f'{wd}/XXXXXX.INP','w',newline='').write(text)
    subprocess.run(['wine','VleachForGUI.exe'],cwd=wd,env=env,capture_output=True,timeout=600)
    for src,dst in (('XXXXXX.INP','input.inp'),('XXXXXX.OUT','output.out'),('XXXXXX.PRF','profiles.prf'),('XXXXXX.PRM','parameters.prm'),('GWIMP.DAT','gwimp.dat')):
        data=open(f'{wd}/{src}','rb').read().replace(b'\r\n',b'\n')
        open(f'{d}/{dst}','wb').write(data)
    print(name, os.path.getsize(f'{d}/output.out'))
# EPA sample
d=f'{OUT}/epa_sample'; shutil.rmtree(d,ignore_errors=True); os.makedirs(d)
S=model / 'Projects' / 'Sample'
for src,dst in (('Sample.inp','input.inp'),('Sample.out','output.out'),('Sample.prf','profiles.prf'),('Sample.prm','parameters.prm'),('Sample_GWIMP.dat','gwimp.dat')):
    open(f'{d}/{dst}','wb').write(open(f'{S}/{src}','rb').read().replace(b'\r\n',b'\n'))
