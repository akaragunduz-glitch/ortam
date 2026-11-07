# CLAUDE.md — Proje çalışma kuralları

Bu dosya, bu depoda çalışan her Claude oturumunun ilk okuduğu belgedir.
Kısa, güncel ve bağlayıcı tutulmalıdır.

## Proje

**site-risk**: Toprak, yeraltı suyu ve toprak gazı için entegre saha
karakterizasyonu, çevresel akıbet ve taşınım modellemesi, insan sağlığı
risk değerlendirmesi ve risk tabanlı ıslah hedefleri yazılımı.

- Amaç ve kapsam: `docs/00-amac-kapsam.md`
- Mimari: `docs/01-mimari.md`
- Modül sözleşmeleri (girdi/çıktı şemaları): `docs/interfaces/`
- Günlük notlar: `docs/gunluk/YYYY-MM-DD.md`
- Kaynak belgeler (EPA, USGS PDF'leri): `docs/kaynaklar/`

## Proje sahibi hakkında

Proje sahibi bilimsel sayısal yöntemler konusunda deneyimlidir ve geçmişte
Fortran ile çalışmıştır; yazılım geliştirme deneyimi sınırlıdır.
- Kodu açıklarken Fortran karşılaştırmaları kullan.
- Her terminal komutunun ne yaptığını tek cümleyle açıkla.
- Yanıtlar Türkçe ve akademik dilde olsun; birimler sembolle yazılsın
  (m3, m3/h, %, mg/kg, µg/L, µg/m3).

## Bağlayıcı kurallar

1. **Denklem hafızadan yazılmaz.** Her denklem `docs/kaynaklar/` altındaki
   bir belgeden alınır ve kodda yorum satırıyla kaynak, sayfa ve denklem
   numarası belirtilir. Örnek: `# RAGS Part B (1991), s. 23, Denk. 4`
2. **Önce test.** Yeni bir hesap fonksiyonundan önce, referans değerlere
   (EPA RSL/VISL hesaplayıcı çıktısı, kılavuz örneği, analitik çözüm)
   dayanan `pytest` testi yazılır.
3. **Birimler açıktır.** Fonksiyon girdileri ve çıktıları birimleriyle
   belgelenir; mümkünse `pint` kullanılır. Değişken adları birim içerebilir
   (`c_soil_mg_kg`).
4. **Tespit altı değerler** sıfır veya DL/2 olarak saklanmaz; ayrı bir
   niteleyici alanla tutulur.
5. **Modüller birbirini doğrudan çağırmaz.** İletişim `docs/interfaces/`
   altındaki Pydantic şemaları ve veritabanı üzerinden olur.
6. **Fortran kodları** `engines/<model>/` altında, sabit derleyici
   sürümüyle Docker içinde derlenir; kaynak kodu değiştirilmişse fark
   (diff) belgelenir.
7. **Plan önce.** Birden fazla dosyayı etkileyen her görevde önce plan
   sunulur, onay beklenir.
8. **Her oturum sonunda** `docs/gunluk/` altına günlük not yazılır:
   yapılanlar, kalanlar, açık sorular.

## Klasör yapısı

```
core/        Bilimsel hesap paketi (arayüzden bağımsız, saf Python)
  exposure/  RAGS maruziyet ve risk denklemleri
  vapor/     Johnson & Ettinger buhar intrüzyonu
  vadose/    Vadoz zon: SSL sızıntı, VLEACH, Richards çözücü (f2py)
  foodchain/ HHRAP bitki ve gıda zinciri
  stats/     ProUCL yöntemleri
  tests/     pytest testleri
engines/     Fortran ve harici kodlar + sarmalayıcılar
  modflow/   MODFLOW 6, MT3D-USGS (FloPy)
  aermod/    AERMOD, AERMET, AERMAP
  vadose/    VLEACH, VS2DT, Richards çözücü Fortran kaynağı
  geochem/   PHREEQC (IPhreeqc), MINTEQA2
api/         FastAPI sunucusu (Faz 2)
worker/      Celery arka plan işçileri (Faz 5)
web/         React + TypeScript arayüzü (Faz 3)
data/reference/  Doğrulama için referans veri setleri
docs/        Belgeler
```

## Komutlar

- Testleri çalıştır: `uv run pytest`
- Bağımlılık ekle: `uv add <paket>`

## Mevcut durum

- Faz 0 başladı (2026-11-07). Klasör yapısı ve kurallar oluşturuldu.
- Sıradaki görev: D0.2 — geliştirme ortamının kurulması ve CI.
