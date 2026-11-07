"""Üç-nokta yöntemiyle hidrolik eğim ve yeraltı suyu akış yönü.

Yöntem: Üç kuyudaki hidrolik yükten geçen düzlem h(x, y) = a·x + b·y + c
çözülür. Hidrolik eğim i = sqrt(a² + b²); akış, eğimin tersi yönündedir.
Kaynak: Klasik üç-nokta problemi (ör. Fetter, Applied Hydrogeology); düzlem
denklemi analitik geometridir. Doğrulama: elle hesap (bkz. testler).

Bu araç MODFLOW gerektirmez; saha verisinden hızlı bir ilk değerlendirme
sağlar. Sayısal akış modeli (MODFLOW 6) aynı modüle ileride eklenecektir.
"""

from __future__ import annotations

import math

import numpy as np
from pydantic import Field

from ortam_core.model import ModelCiktisi, ModelGirdisi, Modul


class Kuyu(ModelGirdisi):
    ad: str
    x_m: float = Field(description="Doğu koordinatı, m (metrik projeksiyon, ör. UTM)")
    y_m: float = Field(description="Kuzey koordinatı, m")
    h_m: float = Field(description="Hidrolik yük, m (kuyu ağzı kotu − suya derinlik)")


class UcNoktaGirdisi(ModelGirdisi):
    kuyular: list[Kuyu] = Field(min_length=3, max_length=3)
    K_m_gun: float | None = Field(None, gt=0, description="Hidrolik iletkenlik, m/gün")
    ne: float | None = Field(None, gt=0, le=1, description="Etkin porozite, -")


class UcNoktaCiktisi(ModelCiktisi):
    egim: float = Field(description="Hidrolik eğim i, m/m")
    akis_azimutu_derece: float = Field(description="Akış yönü, kuzeyden saat yönünde, °")
    darcy_hizi_m_gun: float | None = None
    gozenek_hizi_m_gun: float | None = None


class UcNoktaModulu(Modul):
    ad = "groundwater.uc-nokta"
    baslik = "Yeraltı suyu: üç-nokta akış yönü"
    girdi_tipi = UcNoktaGirdisi
    cikti_tipi = UcNoktaCiktisi

    def calistir(self, girdi: UcNoktaGirdisi) -> UcNoktaCiktisi:
        # Fortran karşılığı: 3x3 doğrusal sistem A·[a, b, c] = h
        A = np.array([[k.x_m, k.y_m, 1.0] for k in girdi.kuyular])
        h = np.array([k.h_m for k in girdi.kuyular])
        if abs(np.linalg.det(A)) < 1e-9:
            raise ValueError("Kuyular aynı doğru üzerinde; üç-nokta çözümü tanımsız.")
        a, b, _ = np.linalg.solve(A, h)

        egim = math.hypot(a, b)
        # Akış -∇h yönünde: doğu bileşeni -a, kuzey bileşeni -b
        azimut = math.degrees(math.atan2(-a, -b)) % 360.0

        uyarilar = []
        if egim < 1e-4:
            uyarilar.append("Eğim çok düşük (< 1e-4); akış yönü belirsiz olabilir.")

        q = v = None
        if girdi.K_m_gun is not None:
            q = girdi.K_m_gun * egim              # Darcy yasası, q = K·i
            if girdi.ne is not None:
                v = q / girdi.ne                  # Ortalama gözenek hızı

        return UcNoktaCiktisi(egim=egim, akis_azimutu_derece=azimut,
                              darcy_hizi_m_gun=q, gozenek_hizi_m_gun=v,
                              uyarilar=uyarilar)

    def ozet(self, c: UcNoktaCiktisi) -> str:
        satirlar = [
            f"Hidrolik eğim       : {c.egim:.5f} m/m",
            f"Akış yönü (azimut)  : {c.akis_azimutu_derece:.1f}° (kuzeyden saat yönünde)",
        ]
        if c.darcy_hizi_m_gun is not None:
            satirlar.append(f"Darcy hızı          : {c.darcy_hizi_m_gun:.4g} m/gün")
        if c.gozenek_hizi_m_gun is not None:
            satirlar.append(f"Gözenek hızı        : {c.gozenek_hizi_m_gun:.4g} m/gün")
        satirlar += [f"UYARI: {u}" for u in c.uyarilar]
        return "\n".join(satirlar)
