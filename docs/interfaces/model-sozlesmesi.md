# Model sözleşmesi

Her ORTAM modülü `ortam_core/model.py` içindeki `Modul` sınıfından türetilir.

| Öğe | Açıklama |
|-----|----------|
| `ad` | Komut satırı adı, `<modül>.<araç>` biçiminde (ör. `groundwater.uc-nokta`, `groundwater.modflow`) |
| `baslik` | İnsan okunur ad |
| `girdi_tipi` | Pydantic sınıfı; her alan birimiyle belgelenir, tanımsız alan hatadır |
| `cikti_tipi` | Pydantic sınıfı; `uyarilar` alanı her çıktıda bulunur |
| `dogrula(veri)` | Ham veriyi tipli girdiye çevirir; hatalı girdi burada reddedilir |
| `calistir(girdi)` | Hesabı yapar; veritabanına ve başka modüllere bağımlı değildir |
| `ozet(cikti)` | Komut satırı için kısa, okunur özet |

## Yeni modül eklerken kontrol listesi

- [ ] `modules/<modül>/<araç>.py` — sınıf, girdi ve çıktı tipleri
- [ ] `modules/<modül>/ornekler/<araç>.yaml` — çalışan örnek girdi
- [ ] `modules/<modül>/tests/test_<araç>.py` — referans testleri
- [ ] `modules/kayit.py` — tek satırlık kayıt
- [ ] `docs/yontem/<modül>.md` — denklemler, parametreler, kaynaklar

Örnek uygulama: `modules/groundwater/uc_nokta.py`.
