# ORTAM — Amaç ve Kapsam

Sürüm: 0.1 · Tarih: 2026-11-07

## 1. Amaç

ORTAM (Open Risk, Transport And Modelling — Ortamlar arası Risk, Taşınım,
Akıbet ve Modelleme), toprak, yeraltı suyu, toprak gazı, hava ve gıda zinciri
ortamları için **modüler bir çevresel modelleme platformudur**.

Platform iki biçimde kullanılır:

1. **Bağımsız modül kullanımı:** Her model tek başına çalıştırılabilir. Örneğin
   kullanıcı yalnızca yeraltı suyu akış ve taşınım modelini, yalnızca AERMOD
   hava dağılımını veya yalnızca jeokimyasal türleşme hesabını kullanabilir.
   Bunun için risk değerlendirmesi, veritabanı veya başka bir modül gerekmez.
2. **Entegre iş akışı kullanımı:** Modüller, Kavramsal Saha Modeli (kaynak →
   taşınım yolu → maruziyet yolu → reseptör) çerçevesinde birbirine bağlanır.
   Örnekleme tasarımından risk tabanlı ıslah hedeflerine uzanan risk zinciri
   bu iş akışlarından biridir; yeni iş akışları (ör. yalnızca yeraltı suyu
   izleme ve trend analizi) eklenebilir.

## 2. Kapsam

### 2.1 Kapsam içi (her modül bağımsız çalışır)

| No | Modül | İçerik |
|----|-------|--------|
| M1 | Proje ve saha | Saha sınırı, koordinat sistemleri, harita altlıkları |
| M2 | Örnekleme tasarımı | Numune sayısı ve konumu (VSP yöntemleri) |
| M3 | Saha ve laboratuvar verisi | Toprak (üst/alt bölme), yeraltı suyu, toprak gazı; GPS koordinatları; QA/QC |
| M4 | Karakterizasyon | Kuyu logu, potansiyometrik yüzey, akış yönü, kriging, 3B jeolojik model |
| M5 | İstatistik | UCL95, arka plan karşılaştırması, tespit altı veri, trend analizi |
| M6 | Yeraltı suyu akışı | MODFLOW 6, MODPATH 7 |
| M7 | Kirletici taşınımı | MODFLOW 6 GWT, MT3D-USGS, Domenico analitik çözümü |
| M8 | Vadoz zon | SSL sızıntı denklemi, VLEACH, Richards çözücü, MODFLOW 6 UZF |
| M9 | Buhar intrüzyonu | Johnson & Ettinger (1991), EPA VISL |
| M10 | Hava dağılımı | AERMOD, AERMET, AERMAP; SSL VF ve PEF |
| M11 | Yüzey akışı ve erozyon | RUSLE, HHRAP denklemleri |
| M12 | Bitki ve gıda zinciri | EPA HHRAP (2005) |
| M13 | Jeokimyasal türleşme | PHREEQC (minteq.v4), MINTEQA2 |
| M14 | Maruziyet ve risk | RAGS A/B/E/F; HQ, HI, ILCR; Monte Carlo |
| M15 | Islah hedefleri ve strateji | Risk tabanlı geri hesap, teknoloji tarama, ıslah hacmi (m3) |
| M16 | Raporlama | PDF ve Word raporları |

### 2.2 Kapsam dışı (sürüm 1.0)

- Ekolojik risk değerlendirmesi (sonraki sürümlerde değerlendirilecek)
- Radyonüklitler
- Çok fazlı akış (NAPL göçünün sayısal modellemesi)
- Gerçek zamanlı sensör verisi akışı

## 3. Hedef kullanıcılar

- Çevre mühendisleri ve hidrojeologlar (saha karakterizasyonu, modelleme)
- Risk değerlendirme uzmanları
- Düzenleyici kurum incelemecileri (salt okuma, rapor inceleme)

## 4. Başarı ölçütleri (sürüm 1.0)

1. Her modül; komut satırından, Python'dan ve web arayüzünden, yalnızca
   kendi girdi dosyasıyla tek başına çalıştırılabilmelidir.
1a. Bir sahanın verisi girildikten sonra risk sonucu ve ıslah hedefi, aynı
   yazılım içinde, kaynağı izlenebilir biçimde üretilebilmelidir.
2. Her hesap modülü, ilgili EPA/USGS referans çıktılarıyla doğrulanmış
   olmalıdır.
3. Her sonuç; kullanılan numuneler, model çalıştırması, parametreler ve
   toksisite tablosu sürümüyle birlikte yeniden üretilebilmelidir.

## 5. Mevzuat ve standart çerçevesi

- Toprak Kirliliğinin Kontrolü ve Noktasal Kaynaklı Kirlenmiş Sahalara Dair
  Yönetmelik ve ilgili teknik rehberler
- US EPA RAGS, RSL, VISL, SSL, HHRAP, ProUCL
- ASTM E1739, ASTM E2081

## 6. Temel varsayımlar ve kısıtlar

- Yazılım Claude ile, modüler olarak geliştirilecektir.
- Zaman kısıtı yoktur; kalite ve doğrulanabilirlik önceliklidir.
- Tüm bileşenler ücretsiz/açık kaynak veya kamu malıdır.
- Mimari web tabanlıdır; yerel (Docker) ve masaüstü (Tauri) çalıştırma
  desteklenir.
