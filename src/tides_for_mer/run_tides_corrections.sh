#!/usr/bin/env bash
set -uo pipefail

# ---------- Configuration ----------
SCRIPT_DIR="/path_script_dir"
PYTHON_SCRIPT="${SCRIPT_DIR}/set_bc_tides_corrections.py"
CONFIG_FILE="${SCRIPT_DIR}/tides_config_new.json"
PYTHON_BIN="/path_python"

# Base path for per-date results folders
RESULTS_BASE="/home/sgaliana/ShareMER_TidalSignal_for_AllDomains/SIC"

# Dates to loop over (YYYYMMDD format)
DATES=("20101224" "20110101" "20201224" "20210101")

# --- Backup the user's hand-edited config BEFORE we touch it ---
BACKUP="${CONFIG_FILE}.domain_backup"
cp "${CONFIG_FILE}" "${BACKUP}"
echo "[OK] Backed up user config to ${BACKUP}"

# Make sure we always restore on exit (even on Ctrl+C or error)
restore() {
	cp "${BACKUP}" "${CONFIG_FILE}"
	rm -f "${BACKUP}"
	echo "[OK] Restored original tides_config_new.json"
}
trap restore EXIT

# ---------- Main loop ----------
for DATE in "${DATES[@]}"; do
	echo "=================================================="
	echo "Processing date: ${DATE}"
	echo "=================================================="
	# Per-date results directory
	RESULTS_DIR="${RESULTS_BASE}/Results_for_MITgcm_SIC_${DATE}"
	mkdir -p "${RESULTS_DIR}"
	
	# Substitute placeholder to date, write the per-date config in place
	sed -i "s/DATEPLACEHOLDER/${DATE}/g" "${CONFIG_FILE}"
	echo "[OK] Substituted DATEPLACEHOLDER to ${DATE} in ${CONFIG_FILE}"

	# Snapshot it into the results folder
	cp "${CONFIG_FILE}" "${RESULTS_DIR}/tides_config_${DATE}.json"
	echo "[OK] Copied config to ${RESULTS_DIR}/tides_config_${DATE}.json"

	# Run the python script
	echo "[RUN] ${PYTHON_BIN} ${PYTHON_SCRIPT}"
        "${PYTHON_BIN}" "${PYTHON_SCRIPT}"	

	# Restore the placeholder version for the NEXT iteration
	cp "${BACKUP}" "${CONFIG_FILE}"

	echo "[DONE] Finished date ${DATE}"
	echo
done

echo "All dates processed successfully."
