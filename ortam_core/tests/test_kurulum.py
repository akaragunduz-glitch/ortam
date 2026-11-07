"""Kurulum ve komut satırı testleri."""

from pathlib import Path

from ortam_core.cli import main


def test_liste(capsys):
    assert main(["liste"]) == 0
    assert "groundwater.uc-nokta" in capsys.readouterr().out


def test_modul_tek_basina_calisir(capsys):
    ornek = Path(__file__).parents[2] / "modules/groundwater/ornekler/uc_nokta.yaml"
    assert main(["calistir", "groundwater.uc-nokta", str(ornek)]) == 0
    assert "Akış yönü" in capsys.readouterr().out
