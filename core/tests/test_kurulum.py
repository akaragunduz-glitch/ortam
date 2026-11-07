"""Kurulum testi: ortamın ve test altyapısının çalıştığını doğrular."""

import numpy as np


def test_numpy_calisiyor():
    # Fortran karşılığı: SUM(A) — NumPy dizileri 0'dan başlar
    a = np.array([1.0, 2.0, 3.0])
    assert a.sum() == 6.0
