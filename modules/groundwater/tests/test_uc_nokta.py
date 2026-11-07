"""Üç-nokta testleri. Referans: elle hesap (düzlem denklemi)."""

import pytest

from modules.groundwater.uc_nokta import UcNoktaModulu


def calistir(kuyular, **kw):
    m = UcNoktaModulu()
    veri = {"kuyular": [dict(zip(("ad", "x_m", "y_m", "h_m"), k)) for k in kuyular], **kw}
    return m.calistir(m.dogrula(veri))


def test_dogu_yonlu_akis():
    # h 100 m'de 1 m düşüyor (doğuya) → i = 0,01, azimut 90°
    c = calistir([("A", 0, 0, 10), ("B", 100, 0, 9), ("C", 0, 100, 10)])
    assert c.egim == pytest.approx(0.01)
    assert c.akis_azimutu_derece == pytest.approx(90.0)


def test_guney_bati_yonlu_akis_ve_hizlar():
    # h kuzeye ve doğuya 0,005 m/m artıyor → akış güneybatıya (225°)
    c = calistir([("A", 0, 0, 10), ("B", 200, 0, 11), ("C", 0, 200, 11)], K_m_gun=10, ne=0.2)
    i = (0.005**2 * 2) ** 0.5
    assert c.egim == pytest.approx(i)
    assert c.akis_azimutu_derece == pytest.approx(225.0)
    assert c.darcy_hizi_m_gun == pytest.approx(10 * i)
    assert c.gozenek_hizi_m_gun == pytest.approx(10 * i / 0.2)


def test_dogrusal_kuyular_hata_verir():
    with pytest.raises(ValueError):
        calistir([("A", 0, 0, 10), ("B", 50, 0, 9), ("C", 100, 0, 8)])


def test_eksik_kuyu_reddedilir():
    with pytest.raises(Exception):
        UcNoktaModulu().dogrula({"kuyular": [{"ad": "A", "x_m": 0, "y_m": 0, "h_m": 1}]})
