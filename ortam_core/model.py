"""ORTAM model arayüzü.

Her modül (yeraltı suyu, vadoz zon, hava, risk ...) bu sözleşmeye uyar.
Böylece her modül:
  1. tek başına (komut satırından, Python'dan veya arayüzden) çalıştırılabilir,
  2. entegre iş akışlarında (ör. risk zinciri) bir halka olarak kullanılabilir.

Fortran karşılaştırması: Bu sınıf, her modelin uyması gereken bir
INTERFACE bloğuna benzer. Her modül aynı üç alt programı sağlar:
dogrula (girdi kontrolü), calistir (hesap), ozet (sonuç raporu).
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from typing import ClassVar

from pydantic import BaseModel


class ModelGirdisi(BaseModel):
    """Tüm modül girdilerinin temel sınıfı. Alanlar birimleriyle belgelenir."""

    model_config = {"extra": "forbid"}  # Tanımsız bir parametre hatadır


class ModelCiktisi(BaseModel):
    """Tüm modül çıktılarının temel sınıfı."""

    uyarilar: list[str] = []


class Modul(ABC):
    """Bağımsız çalışabilen bir ORTAM modülü."""

    ad: ClassVar[str]               # Komut satırı adı, ör. "groundwater"
    baslik: ClassVar[str]           # İnsan okunur ad
    girdi_tipi: ClassVar[type[ModelGirdisi]]
    cikti_tipi: ClassVar[type[ModelCiktisi]]

    def dogrula(self, veri: dict) -> ModelGirdisi:
        """Ham veriyi (YAML, JSON, form) doğrular ve tipli girdiye çevirir."""
        return self.girdi_tipi.model_validate(veri)

    @abstractmethod
    def calistir(self, girdi: ModelGirdisi) -> ModelCiktisi:
        """Modeli çalıştırır. Veritabanına veya başka modüllere bağımlı değildir."""

    def ozet(self, cikti: ModelCiktisi) -> str:
        """Sonucun kısa, okunur özeti (komut satırı için)."""
        return cikti.model_dump_json(indent=2)
