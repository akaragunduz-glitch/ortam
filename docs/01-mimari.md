# Mimari

Sürüm: 0.1 · Tarih: 2026-11-07

Ayrıntılı mimari, teknoloji yığını ve yol haritası için ana belge:
"Entegre Saha Karakterizasyonu, Modelleme ve Risk Değerlendirmesi Yazılımı:
Mimari ve Yol Haritası" (Claude Docs). Bu dosya, kod deposunda Claude'un
başvuracağı özet sürümdür.

## Katmanlar

```
┌──────────────────────────────────────────────────────────────┐
│ Kullanıcı arayüzü (React + TypeScript)                       │
│  2B harita (MapLibre) · 3B (CesiumJS) · grafik/kuyu logu (D3) │
└───────────────────────────┬──────────────────────────────────┘
                            │ HTTPS / JSON
┌───────────────────────────▼──────────────────────────────────┐
│ API sunucusu (FastAPI) — doğrulama, yetki, hızlı hesaplar     │
└──────┬────────────────────┬─────────────────────┬────────────┘
       │                    │                     │
┌──────▼──────┐   ┌─────────▼────────┐   ┌────────▼───────┐
│ PostgreSQL  │   │ İş kuyruğu       │   │ Dosya deposu   │
│ + PostGIS   │   │ (Celery + Redis) │   │ (MinIO)        │
└──────▲──────┘   └─────────┬────────┘   └────────▲───────┘
       │ sonuçlar           │ iş emri             │ dosyalar
┌──────┴────────────────────▼─────────────────────┴───────────┐
│ Hesap motorları (Docker konteynerleri)                        │
│  Python çekirdeği: RAGS, J&E, HHRAP, ProUCL, kriging, vadoz   │
│  Fortran: MODFLOW 6, MT3D-USGS, AERMOD, VLEACH, VS2DT, MINTEQA2│
│  Kütüphaneler: IPhreeqc, GemPy, PyVista                       │
└───────────────────────────────────────────────────────────────┘
```

## Tasarım ilkeleri

1. **Önce bilim, sonra arayüz:** `core/` paketi arayüzden bağımsızdır ve tek
   başına test edilir.
2. **Sözleşme önce:** Her modülün girdi/çıktı şeması `docs/interfaces/`
   altında, kod yazılmadan tanımlanır.
3. **Modüller veritabanı üzerinden konuşur:** Bir modülün çıktısı
   (`model_result`) diğerinin girdisidir.
4. **Fortran kodları dosya tabanlı sarmalayıcıyla çağrılır**; küçük sayısal
   çekirdekler (Richards çözücü) f2py ile bağlanır.
5. **Yeniden üretilebilirlik:** Her model çalıştırmasının girdisi, kod
   sürümü ve parametreleri saklanır.

## Teknoloji yığını (özet)

| Katman | Araç |
|--------|------|
| Dil | Python 3.12, TypeScript |
| Bilimsel | NumPy, SciPy, pandas, statsmodels, PyKrige, GSTools, pint |
| Fortran köprüsü | subprocess, f2py, FloPy |
| API | FastAPI, Pydantic, SQLAlchemy, GeoAlchemy2, Alembic |
| Veritabanı | PostgreSQL + PostGIS |
| İş kuyruğu | Celery + Redis |
| Arayüz | React, MapLibre GL JS, CesiumJS, Plotly, D3.js |
| Paketleme | Docker, Docker Compose |
| Test | pytest, GitHub Actions |
