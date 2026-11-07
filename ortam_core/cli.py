"""ORTAM komut satırı arayüzü.

Kullanım:
  ortam liste                                   Mevcut modülleri listeler
  ortam calistir <modul> <girdi.yaml>           Bir modülü tek başına çalıştırır
  ortam calistir <modul> <girdi.yaml> --json    Sonucu JSON olarak yazar
"""

import argparse
import sys
from pathlib import Path

import yaml

from modules.kayit import MODULLER, modul_al


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(prog="ortam", description="ORTAM çevresel modelleme platformu")
    alt = p.add_subparsers(dest="komut", required=True)
    alt.add_parser("liste", help="Mevcut modülleri listele")
    c = alt.add_parser("calistir", help="Bir modülü tek başına çalıştır")
    c.add_argument("modul")
    c.add_argument("girdi", type=Path)
    c.add_argument("--json", action="store_true", help="Sonucu JSON olarak yaz")
    a = p.parse_args(argv)

    if a.komut == "liste":
        for ad in MODULLER:
            print(f"{ad:28s} {modul_al(ad).baslik}")
        return 0

    modul = modul_al(a.modul)
    girdi = modul.dogrula(yaml.safe_load(a.girdi.read_text(encoding="utf-8")))
    cikti = modul.calistir(girdi)
    print(cikti.model_dump_json(indent=2) if a.json else modul.ozet(cikti))
    return 0


if __name__ == "__main__":
    sys.exit(main())
