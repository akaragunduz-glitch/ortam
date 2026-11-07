# CLAUDE.md — Proje çalışma kuralları

Bu dosya, bu depoda çalışan her Claude oturumunun ilk okuduğu belgedir.
Kısa, güncel ve bağlayıcı tutulmalıdır.

## Proje

**ORTAM** (Open Risk, Transport And Modelling — Ortamlar arası Risk,
Taşınım, Akıbet ve Modelleme): modüler bir çevresel modelleme platformu.
Yeraltı suyu, vadoz zon, buhar intrüzyonu, hava dağılımı, yüzey akışı, gıda
zinciri, jeokimya, istatistik, karakterizasyon, risk ve ıslah modülleri
içerir.

**Temel ilke: her modül tek başına çalışır.** Kullanıcı yalnızca yeraltı suyu
modelini kullanabilir; risk değerlendirmesi, modülleri birbirine bağlayan
iş akışlarından yalnızca biridir (`workflows/risk_zinciri/`).

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
5. **Modüller birbirini doğrudan çağırmaz.** Her modül `ortam_core.model.Modul`
   sözleşmesine uyar (`dogrula`, `calistir`, `ozet`) ve veritabanı olmadan,
   yalnızca bir girdi dosyasıyla çalışabilmelidir. Modüller arası bağlantı
   yalnızca `workflows/` altındaki iş akışlarında kurulur.
5a. **Her modül üç yoldan çalıştırılabilir olmalıdır:** komut satırı
   (`ortam calistir <modul> <girdi.yaml>`), Python (`modul.calistir(...)`)
   ve web arayüzü. Her modülde `ornekler/` altında en az bir örnek girdi
   ve `tests/` altında referans testleri bulunur.
6. **Fortran kodları** `engines/<model>/` altında, sabit derleyici
   sürümüyle Docker içinde derlenir; kaynak kodu değiştirilmişse fark
   (diff) belgelenir.
7. **Plan önce.** Birden fazla dosyayı etkileyen her görevde önce plan
   sunulur, onay beklenir.
8. **Her oturum sonunda** `docs/gunluk/` altına günlük not yazılır:
   yapılanlar, kalanlar, açık sorular.

## Klasör yapısı

```
ortam_core/        Ortak çekirdek: model sözleşmesi, birimler, CLI
modules/           Bağımsız modüller (her biri tek başına çalışır)
  groundwater/     Yeraltı suyu: üç-nokta, MODFLOW 6 akış ve GWT taşınım
  vadose/          Vadoz zon: SSL sızıntı, VLEACH, Richards çözücü, UZF
  vapor/           Buhar intrüzyonu: Johnson & Ettinger
  air/             Hava dağılımı: AERMOD, emisyon faktörleri
  surface/         Yüzey akışı ve erozyon: RUSLE
  foodchain/       Bitki ve gıda zinciri: HHRAP
  geochem/         Jeokimya: PHREEQC, MINTEQA2
  stats/           İstatistik: ProUCL yöntemleri
  characterization/ Kuyu logu, potansiyometrik yüzey, kriging
  sampling/        Örnekleme tasarımı (VSP yöntemleri)
  risk/            Maruziyet ve risk (RAGS)
  remediation/     Risk tabanlı ıslah hedefleri ve strateji
  kayit.py         Modül kaydı (CLI ve arayüz buradan okur)
  <modul>/ornekler/  Örnek girdi dosyaları (YAML)
  <modul>/tests/     Referans testleri
workflows/         Modülleri birbirine bağlayan iş akışları
  risk_zinciri/    Örnekleme → modelleme → risk → ıslah
engines/           Fortran ve harici kodlar (MODFLOW, AERMOD, VLEACH, ...)
api/               FastAPI sunucusu (Faz 2)
worker/            Celery arka plan işçileri (Faz 5)
web/               React + TypeScript arayüzü (Faz 3)
data/reference/    Doğrulama için referans veri setleri
docs/              Belgeler
```

## Komutlar

- Testleri çalıştır: `uv run pytest`
- Modülleri listele: `uv run ortam liste`
- Bir modülü tek başına çalıştır:
  `uv run ortam calistir groundwater.uc-nokta modules/groundwater/ornekler/uc_nokta.yaml`
- Bağımlılık ekle: `uv add <paket>`

## Mevcut durum

- Faz 0 başladı (2026-11-07). Klasör yapısı, kurallar, model sözleşmesi ve
  CLI oluşturuldu. İlk bağımsız araç: `groundwater.uc-nokta` (testli).
- Sıradaki görev: D0.2 — geliştirme ortamının kurulması ve CI.
