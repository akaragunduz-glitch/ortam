# ORTAM

**Open Risk, Transport And Modelling** — Ortamlar arası Risk, Taşınım, Akıbet
ve Modelleme platformu.

Toprak, yeraltı suyu, toprak gazı, hava ve gıda zinciri için modüler çevresel
modelleme yazılımı. Her modül (yeraltı suyu, vadoz zon, buhar intrüzyonu, hava
dağılımı, jeokimya, istatistik, risk, ıslah ...) **tek başına** çalıştırılabilir
veya iş akışlarıyla birbirine bağlanabilir.

```bash
uv run ortam liste
uv run ortam calistir groundwater.uc-nokta modules/groundwater/ornekler/uc_nokta.yaml
```

- Amaç ve kapsam: [docs/00-amac-kapsam.md](docs/00-amac-kapsam.md)
- Mimari: [docs/01-mimari.md](docs/01-mimari.md)
- Model sözleşmesi: [docs/interfaces/model-sozlesmesi.md](docs/interfaces/model-sozlesmesi.md)
- Claude çalışma kuralları: [CLAUDE.md](CLAUDE.md)

Durum: Faz 0 — altyapı kurulumu.
