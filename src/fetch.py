"""Download FAERS adverse event reports (openFDA) and keep, per report, the list of drugs, patient age, and whether the patient died."""
import requests, pandas as pd, time, os
S, E = int(os.environ.get('S', 0)), int(os.environ.get('E', 15000))   # 100 reports per call, run in chunks
rows = []
for skip in range(S, E, 100):
    r = requests.get("https://api.fda.gov/drug/event.json", params={
        "search": "receivedate:[20230101 TO 20231231] AND patient.patientonsetage:[18 TO 120] AND _exists_:patient.drug.openfda.generic_name",
        "limit": 100, "skip": skip}, timeout=120)
    if r.status_code != 200: print("stop", skip, r.status_code); break
    for x in r.json()["results"]:
        p = x.get("patient", {}); names = []
        for dr in p.get("drug", []):
            g = dr.get("openfda", {}).get("generic_name")
            if g: names.append(g[0].lower())
        rows.append({"id": x["safetyreportid"], "age": p.get("patientonsetage"), "sex": p.get("patientsex"),
                     "died": int(x.get("seriousnessdeath") == "1"), "n_drugs_listed": len(p.get("drug", [])),
                     "drugs": "|".join(sorted(set(names)))})
    time.sleep(0.2)
out = pd.DataFrame(rows)
if os.path.exists("data/reports.csv") and os.path.getsize("data/reports.csv") > 10 and S > 0: out = pd.concat([pd.read_csv("data/reports.csv"), out])
out.to_csv("data/reports.csv", index=False); print(len(out))
