#!/usr/bin/env bash
set -euo pipefail

# Prepare a paper plasma-expansion case:
#   distribution  : Maxwellian, Kappa, or polytropic electrons
#   species       : electron + one ion species
#   mass ratio    : configurable mi/me, Z = 1
#   domain        : [0, 1024 lambda_D0] x [0, 4 lambda_D0]
#   initial slab  : [0, 128 lambda_D0] x [0, 4 lambda_D0]
#   mesh          : dx = dy = 1.0 or 0.5 lambda_D0
#   time step     : dt = 0.05 or 0.025 / omega_pe
#
# Usage:
#   bash scripts/setup_maxwell_mi400_case.sh [ppc] [nt] [boundary] [seed] [dt] [dx] [distribution] [shape] [mi/me]
#
# Examples:
#   bash scripts/setup_maxwell_mi400_case.sh 1000 20000 thermal 101 0.05 1.0
#   bash scripts/setup_maxwell_mi400_case.sh 1000 20000 thermal 101 0.05 1.0 kappa 2 400
#   bash scripts/setup_maxwell_mi400_case.sh 1000 20000 thermal 101 0.05 1.0 polytropic 2 400
#
# The paper-level ppc=80000 case is large. Use a short ppc=10, nt=200 run
# first. Runs shorter than omega_pi*t=50 are treated as smoke tests and are
# deliberately not archived as paper results.

ppg="${1:-1000}"
nt="${2:-20000}"
left_boundary="${3:-thermal}"
random_seed="${4:-101}"
dt="${5:-0.05}"
dx="${6:-1.0}"
distribution="${7:-maxwellian}"
distribution_shape="${8:-}"
mass_ratio="${9:-400}"

distribution="$(printf '%s' "$distribution" | tr '[:upper:]' '[:lower:]')"
case "$distribution" in
  maxwellian|maxwell)
    distribution="maxwellian"
    distribution_tag="maxwellian"
    distribution_shape="none"
    electron_kappa="2.0"
    polytropic_gamma="2.0"
    ;;
  kappa)
    distribution_shape="${distribution_shape:-2}"
    if ! [[ "$distribution_shape" =~ ^[0-9]+([.][0-9]+)?$ ]] || \
       ! awk -v value="$distribution_shape" 'BEGIN { exit !(value > 1.5) }'; then
      echo "Kappa shape must be a number greater than 1.5" >&2
      exit 2
    fi
    shape_tag="${distribution_shape//./p}"
    distribution_tag="kappa${shape_tag}"
    electron_kappa="$distribution_shape"
    polytropic_gamma="2.0"
    ;;
  polytropic|poly)
    distribution="polytropic"
    distribution_shape="${distribution_shape:-2}"
    if ! [[ "$distribution_shape" =~ ^[0-9]+([.][0-9]+)?$ ]] || \
       ! awk -v value="$distribution_shape" 'BEGIN { exit !(value > 1.0 && value <= 3.0) }'; then
      echo "Polytropic gamma must satisfy 1 < gamma <= 3" >&2
      exit 2
    fi
    shape_tag="${distribution_shape//./p}"
    distribution_tag="poly${shape_tag}"
    electron_kappa="2.0"
    polytropic_gamma="$distribution_shape"
    ;;
  *)
    echo "distribution must be 'maxwellian', 'kappa', or 'polytropic'" >&2
    exit 2
    ;;
esac

case "$left_boundary" in
  thermal|specular) ;;
  *)
    echo "left boundary must be 'thermal' or 'specular'" >&2
    exit 2
    ;;
esac

case "$dt" in
  0.05|.05) dt="0.05"; dt_tag="005" ;;
  0.025|.025) dt="0.025"; dt_tag="0025" ;;
  *)
    echo "dt must be 0.05 or 0.025" >&2
    exit 2
    ;;
esac

case "$dx" in
  1|1.0) dx="1.0"; dx_tag="1"; mesh_nx=1025; mesh_ny=5; slab_cells_x=128; slab_cells_y=4 ;;
  0.5|.5) dx="0.5"; dx_tag="05"; mesh_nx=2049; mesh_ny=9; slab_cells_x=256; slab_cells_y=8 ;;
  *)
    echo "dx must be 1.0 or 0.5" >&2
    exit 2
    ;;
esac

if ! [[ "$ppg" =~ ^[1-9][0-9]*$ && "$nt" =~ ^[1-9][0-9]*$ && \
        "$random_seed" =~ ^[0-9]+$ && "$mass_ratio" =~ ^[1-9][0-9]*$ ]]; then
  echo "ppc, nt, and mi/me must be positive integers; seed must be a non-negative integer" >&2
  exit 2
fi

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
app_dir="$repo_root/PIC-IFE_GEC"
archive_root="${PIC_ARCHIVE_ROOT:-$repo_root/verification_runs}"
mcc_file="$app_dir/MCC_jw/code/Interface_IFE/MCCInterface.f90"

species_count=2
initial_particles=$((ppg * slab_cells_x * slab_cells_y * species_count))
particle_capacity=$(((initial_particles * 12 + 9) / 10))
if (( particle_capacity < 100000 )); then
  particle_capacity=100000
fi

# Resource model calibrated from the completed dx=1, dt=0.05, ppc=40000 run:
# 40960000 initial particles, nt=20000, 75036 s elapsed, 5522 MiB MaxRSS.
# Add 35% memory and 50% time headroom, then round to scheduler-friendly values.
reference_particles=40960000
reference_nt=20000
reference_elapsed_seconds=75036
reference_maxrss_mib=5522

estimated_memory_mib=$(( \
  (reference_maxrss_mib * initial_particles + reference_particles - 1) \
  / reference_particles \
))
requested_memory_mib=$(((estimated_memory_mib * 135 + 99) / 100))
requested_memory_mib=$((((requested_memory_mib + 499) / 500) * 500))
if (( requested_memory_mib < 1000 )); then
  requested_memory_mib=1000
fi

estimated_runtime_seconds=$(( \
  (reference_elapsed_seconds * initial_particles * nt \
    + reference_particles * reference_nt - 1) \
  / (reference_particles * reference_nt) \
))
requested_runtime_seconds=$(((estimated_runtime_seconds * 3 + 1) / 2))
requested_walltime_hours=$(((requested_runtime_seconds + 3599) / 3600))
if (( requested_walltime_hours < 12 )); then
  requested_walltime_hours=12
fi
requested_walltime_hours=$((((requested_walltime_hours + 11) / 12) * 12))
requested_walltime_days=$((requested_walltime_hours / 24))
requested_walltime_remainder=$((requested_walltime_hours % 24))
if (( requested_walltime_days > 0 )); then
  printf -v requested_walltime '%d-%02d:00:00' \
    "$requested_walltime_days" "$requested_walltime_remainder"
else
  printf -v requested_walltime '%02d:00:00' "$requested_walltime_hours"
fi

target_t30_step="$(awk -v dt="$dt" -v mi="$mass_ratio" 'BEGIN { printf "%.0f", 30.0 * sqrt(mi) / dt }')"
target_t50_step="$(awk -v dt="$dt" -v mi="$mass_ratio" 'BEGIN { printf "%.0f", 50.0 * sqrt(mi) / dt }')"
printf -v target_t30_file_step '%06d' "$target_t30_step"
printf -v target_t50_file_step '%06d' "$target_t50_step"
if (( nt >= target_t50_step )); then
  run_mode="production"
  archive_results="yes"
else
  run_mode="smoke"
  archive_results="no"
fi

backup_once() {
  local f="$1"
  if [[ -f "$f" && ! -f "$f.maxwell_mi400.bak" ]]; then
    cp "$f" "$f.maxwell_mi400.bak"
  fi
}

replace_block() {
  local file="$1"
  local perl_expr="$2"
  local tmp
  tmp="$(mktemp)"
  perl -0pe "$perl_expr" "$file" > "$tmp"
  mv "$tmp" "$file"
}

backup_once "$app_dir/code/PIC/Main_IFE_Test_2.f90"
backup_once "$mcc_file"
backup_once "$app_dir/code/Data/PIC_MAIN_PARAM_2D.f90"
backup_once "$app_dir/code/In-Output/Output_To_Tecplot_IJK_2D.f90"
backup_once "$app_dir/OUTPUT_velocity.f90"
backup_once "$app_dir/Output_Energy.f90"
backup_once "$app_dir/INPUT/mesh.inp"
backup_once "$app_dir/INPUT/object.inp"
backup_once "$app_dir/INPUT/pic.inp"
backup_once "$app_dir/MCC_jw/input/gas.txt"
backup_once "$app_dir/MCC_jw/input/controlflow.txt"

# Initial slab [0,128] x [0,4], embedded in full [0,1024] x [0,4] domain.
perl -0pi -e 's/dxmaxmax = [0-9.]+/dxmaxmax = 128.0/' "$app_dir/code/PIC/Main_IFE_Test_2.f90"
perl -0pi -e 's/dymaxmax = [0-9.]+/dymaxmax = ymax/' "$app_dir/code/PIC/Main_IFE_Test_2.f90"
perl -0pi -e 's/    call OUTPUT_velocity\(it\)\n    call Output_Energy\(it\)/    If (Mod(it,1000).eq.0 .Or. it == 1 .Or. it == nt) Then\n        call OUTPUT_velocity(it)\n        call Output_Energy(it)\n    Endif/' "$app_dir/code/PIC/Main_IFE_Test_2.f90"
perl -0pi -e 's/    If \(Mod\(it,1000\)\.eq\.0 \.Or\. it == 1\) Then/    If (Mod(it,1000).eq.0 .Or. it == 1 .Or. it == nt) Then/' "$app_dir/code/PIC/Main_IFE_Test_2.f90"

# Keep the global particle sanity check consistent with the selected ppc.
perl -0pi -e "s/N_part_max = [0-9]+/N_part_max = $particle_capacity/" "$app_dir/code/Data/PIC_MAIN_PARAM_2D.f90"

# Force the MCC/JW particle bundle to use species data read from INPUT/pic.inp.
perl -0pi -e 's/    Use Field_2D, Only:dens0/    Use Field_2D, Only: dens0/' "$mcc_file"
if ! grep -q 'PIC_EXPANSION_MI400' "$mcc_file"; then
  perl -0pi -e 's/(        Call GasInitPegasus\(ControlFlowGlobal\).*\n)/$1        ! PIC_EXPANSION_MI400: use the two species read from INPUT\/pic.inp.\n        ControlFlowGlobal%Ns = Min(ControlFlowGlobal%Ns, ispe_tot - 1)\n        Do i = 0, ControlFlowGlobal%Ns\n            SpecyGlobal(i)%Charge = qs(i+1) * q_ref\n            SpecyGlobal(i)%Mass = xm(i+1) * m_ref\n            If (Allocated(tmpj)) Then\n                SpecyGlobal(i)%InitTemperature = tmpj(i+1) * T_ref\n                SpecyGlobal(i)%Temperature = SpecyGlobal(i)%InitTemperature\n            End If\n            If (Allocated(dens0)) Then\n                SpecyGlobal(i)%InitDensity = dens0(i+1)\n                SpecyGlobal(i)%Density = SpecyGlobal(i)%InitDensity\n            End If\n        End Do\n/' "$mcc_file"
fi

# Let ParticlePerGrid in INPUT/pic.inp control the initial slab particle count.
perl -0pi -e 's/        !======= for wsy paper case ==========\n        PB%NParNormal = 6000000\n        NPArMax = PB%NParNormal\n        If\(Allocated\(PB%PO\)\) Deallocate\(PB%PO\)\n        Allocate\(PB%PO\(NPArMax\)\)\n        !======================================\n/        ! Initial particle storage is allocated below from ParticlePerGrid.\n/' "$mcc_file"
perl -0pi -e 's/           PB%NPar = 5000 \* \(dxmaxmax-dxminmin\) \* \(dymaxmax-dyminmin\)\s*\n            PB%Weight = affp_bjw\(isp\)\s*\n            !PB%Weight = dens0\(isp\)\*n_ref \* RegionVolume \/ PB%NPar\s*\n            !print\*,dens0\(isp\)\*n_ref\s*\n            PB%NParNormal = 5000 \* \(dxmaxmax-dxminmin\) \* \(dymaxmax-dyminmin\)[^\n]*\n/            If (delta_global == 0) Then\n                RegionVolume = (dxmaxmax-dxminmin)*L_ref * (dymaxmax-dyminmin)*L_ref\n            Elseif (delta_global == 1) Then\n                RegionVolume = (dxmaxmax-dxminmin)*L_ref * PI*(dymaxmax**2-dyminmin**2)*L_ref**2\n            Endif\n            PB%NPar = INT(DBLE(ParticlePerGrid) * (dxmaxmax-dxminmin) * (dymaxmax-dyminmin))\n            PB%Weight = dens0(isp)*n_ref * RegionVolume \/ DBLE(PB%NPar)\n            PB%NParNormal = PB%NPar\n/s' "$mcc_file"
perl -0pi -e 's/NPArMax = Ceiling\(3\.0 \* PB%NParNormal\)/NPArMax = Ceiling(1.2D0 * PB%NParNormal)/' "$mcc_file"
perl -0pi -e 's/PB%NPar = INT\(DBLE\(ParticlePerGrid\) \* \(dxmaxmax-dxminmin\) \* \(dymaxmax-dyminmin\)\)/PB%NPar = INT(DBLE(ParticlePerGrid) * (dxmaxmax-dxminmin) * (dymaxmax-dyminmin) \/ (hx(1)*hx(2)))/' "$mcc_file"

# y boundaries: periodic particle wrap, matching the 1D-dominant slab setup.
y_block='        Elseif (PO%Y < dymin) Then
            If (delta == 0) Then
                Do While (PO%Y < dymin)
                    PO%Y = PO%Y + (dymax - dymin)
                End Do
            Elseif (delta == 1) Then
                Print*, "wrong cross"
                Print*, PO
                pause
            Endif
        Elseif (PO%Y > dymax) Then
            If (delta == 0) Then
                Do While (PO%Y > dymax)
                    PO%Y = PO%Y - (dymax - dymin)
                End Do
            Endif
        Endif'
export PAPER_Y_BOUNDARY_BLOCK="$y_block"
replace_block "$mcc_file" 'BEGIN { $r = $ENV{"PAPER_Y_BOUNDARY_BLOCK"} . "\n" } s/        Elseif \(PO%Y < dymin\) Then.*?        If \(N_objects/$r        If (N_objects/s'
unset PAPER_Y_BOUNDARY_BLOCK

# Two-species output: avoid historical third-species columns.
perl -0pi -e 's/WRITE\(1,\('\''\(A150\)'\''\)\) '\''VARIABLES = "x" "y" "R" "Phi" "Rho" "Rho_ele" "Rho_H" "Rho_C" "efx" "efy" "Ek_ele" "Ek_ion" "Ek_tot_ele" "Ek_tot_ion" "node_type1"  "node_type2"'\''/WRITE(1,('\''(A150)'\'')) '\''VARIABLES = "x" "y" "R" "Phi" "Rho" "Rho_ele" "Rho_ion" "efx" "efy" "Ek_ele" "Ek_ion" "Ek_tot_ele" "Ek_tot_ion" "node_type1"  "node_type2"'\''/' "$app_dir/code/In-Output/Output_To_Tecplot_IJK_2D.f90"
perl -0pi -e 's/          WRITE\(1,50\) HP\(1:2,i\), radius\(i\),Phi\(i,1\), Rho\(i,1\), Rho_s\(i,1,1\), Rho_s\(i,1,2\),Rho_s\(i,1,3\), efx\(i, 1\), efy\(i, 1\), Ek_s\(i,1,1\), Ek_s\(i,1,2\), Ek_tot\(i,1,1\), Ek_tot\(i,1,2\), node_type\(1, i\), node_type\(2, i\)/          WRITE(1,*) HP(1:2,i), radius(i), Phi(i,1), Rho(i,1), Rho_s(i,1,1), Rho_s(i,1,2), \&\n                     efx(i,1), efy(i,1), Ek_s(i,1,1), Ek_s(i,1,2), Ek_tot(i,1,1), Ek_tot(i,1,2), \&\n                     node_type(1,i), node_type(2,i)/' "$app_dir/code/In-Output/Output_To_Tecplot_IJK_2D.f90"

perl -0pi -e 's/        do num=1,N\n            !WRITE\(544, '\''\(A150\)'\''\) [^\n]*\n            WRITE\(544, 18\) driftVelocity\(1,num, 1\), driftVelocity\(2,num, 1\), driftVelocity\(3,num, 1\), thermalVelocity\(1,num, 1\), thermalVelocity\(2,num, 1\), thermalVelocity\(3,num, 1\),x_num\(1,num\),x_num\(2,num\),x_num\(3,num\),num\n        end do\n    \n        \n18      FORMAT[^\n]*\n/        WRITE(544, '\''(A)'\'') '\''# drift_e drift_i thermal_e thermal_i count_e count_i cell'\''\n        do num=1,N\n            WRITE(544, *) driftVelocity(1,num,1), driftVelocity(2,num,1), thermalVelocity(1,num,1), \&\n                          thermalVelocity(2,num,1), x_num(1,num), x_num(2,num), num\n        end do\n    \n        \n/s' "$app_dir/OUTPUT_velocity.f90"
perl -0pi -e 's/WRITE\(544, \*\) driftVelocity/WRITE(544, '\''(6(ES16.8,1X),I8)'\'') driftVelocity/' "$app_dir/OUTPUT_velocity.f90"
perl -0pi -e 's/    If \(Mod\(it,1000\)\.eq\.0 \.Or\. it == 1\) Then/    If (Mod(it,1000).eq.0 .Or. it == 1 .Or. it == nt) Then/' "$app_dir/OUTPUT_velocity.f90"
if ! grep -q 'Use TimeControl, Only: nt' "$app_dir/OUTPUT_velocity.f90"; then
  perl -0pi -e 's/(Use Constant_Variable_2D\r?\n)/$1Use TimeControl, Only: nt\n/' "$app_dir/OUTPUT_velocity.f90"
fi

perl -0pi -e 's/        WRITE\(543, 17\) \&\n            dN_dE\(1, i\), dN_dE\(2, i\),dN_dE\(3, i\), Ek_per\(1, i\), Vk_per\(1, i\), Vk_per\(2, i\), Vk_per\(3, i\), dN_dV\(1, i\), dN_dV\(2, i\),dN_dV\(3, i\)/        WRITE(543, *) dN_dE(1,i), dN_dE(2,i), Ek_per(1,i), Vk_per(1,i), Vk_per(2,i), dN_dV(1,i), dN_dV(2,i)/' "$app_dir/Output_Energy.f90"
perl -0pi -e 's/\n    17 FORMAT[^\n]*\n/\n/' "$app_dir/Output_Energy.f90"
perl -0pi -e 's/WRITE\(543,\s*17\)/WRITE(543, *)/g; s/\r?\n\s*17 FORMAT[^\r\n]*(\r?\n)/$1/g' "$app_dir/Output_Energy.f90"

cat > "$app_dir/INPUT/mesh.inp" <<'EOF_MESH'
0., 0., 0.
1024.0, 4.0, 0.
EOF_MESH
cat >> "$app_dir/INPUT/mesh.inp" <<EOF_MESH
$mesh_nx, $mesh_ny, 0
$dx, $dx
1.0
EOF_MESH

cat > "$app_dir/INPUT/object.inp" <<'EOF_OBJECT'
0,-2
4
0, 0.0, 1024.0, 0.0, 1024.0, 4.0
0, 0.0,    0.0, 0.0, 1024.0, 0.0
0, 0.0,    0.0, 0.0,    0.0, 4.0
0, 0.0,    0.0, 4.0, 1024.0, 4.0
EOF_OBJECT

ion_mass="$(awk -v ratio="$mass_ratio" 'BEGIN { printf "%.10E", 9.1095e-31 * ratio }')"
cat > "$app_dir/INPUT/pic.inp" <<EOF_PIC
! ${distribution_tag} mi/me=${mass_ratio} two-species plasma expansion
0, 0
.false.
.false.
1.0D21
1.0
$ppg
1.0D0, 1.0D0, 1.0D0
$particle_capacity
.FALSE., .FALSE.
2
$nt, $dt
1000000
1
1000000
1, 1
2000
1000
5000
1, 1, 1
1.0, 1.0, 1.0
4.0
-1.6022D-19, 9.1095D-31
 1.6022D-19, $ion_mass
.FALSE., .FALSE., .FALSE., .FALSE.
.TRUE., .TRUE., .TRUE., .TRUE.
.FALSE., .FALSE., .TRUE., .TRUE., .FALSE., .FALSE.
.FALSE., .TRUE., .FALSE., .FALSE., .FALSE., .FALSE.
.TRUE., .FALSE., .FALSE., .FALSE., .FALSE., .FALSE.
.FALSE., .FALSE., .FALSE., .FALSE., .FALSE., .FALSE.
0.0, 0.0, 0.0, 0.0
1.0, 1.0, 0.0
2
1
1.0D21
0.0
1.0
0.0, 0.0, 0.0
0.0, 0.0, 0.0
128.0
0.0
0.0
2
1.0D21
0.0
0.01
0.0, 0.0, 0.0
0.0, 0.0, 0.0
128.0
0.0
0.0
3.0
EOF_PIC

cat > "$app_dir/MCC_jw/input/gas.txt" <<'EOF_GAS'
&Ngfile
Ng           =   1,
GasName      =   "He",
PGas         =   0.0,
TGas         =   11605.0,
TElectron    =   11605.0,
GasRatio     =   1.0,
/
EOF_GAS

cat > "$app_dir/MCC_jw/input/controlflow.txt" <<EOF_CF
&ControlFlow
CF%InitDensity        =   1.0d21,
CF%NRun               =   $nt,
CF%NDiagShort         =   1000,
CF%NDiagLong          =   5000,
CF%ParticlePerGrid    =   $ppg,
CF%withRecombination  =   .false.,
/
EOF_CF

mkdir -p "$app_dir/OUTPUT/Field" "$app_dir/OUTPUT/Velocity" "$app_dir/OUTPUT/Particle" \
         "$app_dir/OUTPUT/Global" "$app_dir/OUTPUT/Phase" "$app_dir/OUTPUT/Energy" \
         "$app_dir/OUTPUT/History" "$app_dir/OUTPUT/Average" "$app_dir/DUMP"

case_name="${distribution_tag}_mi${mass_ratio}_dx${dx_tag}_dt${dt_tag}_ppc${ppg}_seed${random_seed}_${left_boundary}"

cat > "$app_dir/case_config.txt" <<EOF_CONFIG
case_name = ${case_name}
distribution = ${distribution}
distribution_parameter = ${distribution_shape}
electron_kappa = ${electron_kappa}
electron_polytropic_gamma = ${polytropic_gamma}
left_boundary = ${left_boundary}
species = electron, ion
ion_mass_over_electron_mass = ${mass_ratio}
ion_charge_state = 1
density_ref_m3 = 1.0e21
Te0_eV = 1.0
Ti0_eV = 0.01
domain_lambdaD = [0,1024] x [0,4]
initial_slab_lambdaD = [0,128] x [0,4]
mesh = ${mesh_nx} x ${mesh_ny}
dx_lambdaD = ${dx}
dy_lambdaD = ${dx}
dt_omega_pe = ${dt}
nt = ${nt}
run_mode = ${run_mode}
archive_results = ${archive_results}
particles_per_cell_per_species = ${ppg}
initial_particles_total = ${initial_particles}
particle_capacity = ${particle_capacity}
slurm_cpus = 1
slurm_memory_mib = ${requested_memory_mib}
slurm_walltime = ${requested_walltime}
target_omega_pi_t_30_step = ${target_t30_step}
target_omega_pi_t_50_step = ${target_t50_step}
diagnostic_stride_global = 1000
field_velocity_output_stride = 1000
random_seed = ${random_seed}
archive_root = ${archive_root}
EOF_CONFIG

cat > "$app_dir/run_${case_name}.slurm" <<EOF_SLURM
#!/bin/bash
#SBATCH -J ${distribution_tag}_${mass_ratio}_${ppg}
#SBATCH -p comp
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --mem=${requested_memory_mib}M
#SBATCH -t ${requested_walltime}
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.err

module purge
module load intel/2022.1
module load cmake/3.23.5

cd "$app_dir"
ulimit -s unlimited
export PIC_RANDOM_SEED=${random_seed}
export PIC_ELECTRON_DISTRIBUTION=${distribution}
export PIC_KAPPA=${electron_kappa}
export PIC_POLYTROPIC_GAMMA=${polytropic_gamma}
export PIC_LEFT_BOUNDARY=${left_boundary}

if command -v flock >/dev/null 2>&1; then
  exec 9>.pic_run.lock
  if ! flock -n 9; then
    echo "another PIC case is already using this working directory" >&2
    exit 3
  fi
fi

rm -f run.log OUTPUT/global_diagnostics.csv
mkdir -p OUTPUT/{Field,Velocity,Particle,Global,Phase,Energy,History,Average} DUMP

run_start_utc="\$(date -u +%Y-%m-%dT%H:%M:%SZ)"
run_start_epoch="\$(date +%s)"
set +e
./1DPIC > run.log 2>&1
run_status="\$?"
set -e
run_end_utc="\$(date -u +%Y-%m-%dT%H:%M:%SZ)"
run_end_epoch="\$(date +%s)"
run_elapsed_seconds="\$((run_end_epoch - run_start_epoch))"
git_revision="\$(git -C "$repo_root" rev-parse HEAD 2>/dev/null || echo unknown)"

cat > run_metadata.txt <<EOF_METADATA
case_name = ${case_name}
distribution = ${distribution}
distribution_parameter = ${distribution_shape}
left_boundary = ${left_boundary}
ion_mass_over_electron_mass = ${mass_ratio}
run_mode = ${run_mode}
git_revision = \$git_revision
slurm_job_id = \${SLURM_JOB_ID:-unknown}
hostname = \$(hostname)
run_start_utc = \$run_start_utc
run_end_utc = \$run_end_utc
run_elapsed_seconds = \$run_elapsed_seconds
exit_status = \$run_status
EOF_METADATA

if [[ "\$run_status" -ne 0 ]]; then
  echo "1DPIC failed with exit status \$run_status" >&2
  exit "\$run_status"
fi

if [[ "${archive_results}" == "yes" ]]; then
  PIC_APP_DIR="$app_dir" PIC_ARCHIVE_ROOT="$archive_root" \
    bash "$repo_root/scripts/archive_verification_case.sh"
else
  echo "smoke test completed; results remain in $app_dir/OUTPUT and are not archived"
fi
EOF_SLURM

pic_ispe="$(awk 'NR==11 {print $1}' "$app_dir/INPUT/pic.inp" | tr -d ',')"
pic_nt="$(awk 'NR==12 {print $1}' "$app_dir/INPUT/pic.inp" | tr -d ',')"
if [[ "$pic_ispe" != "2" || "$pic_nt" != "$nt" ]]; then
  echo "case validation failed: ispe_tot=$pic_ispe nt=$pic_nt" >&2
  exit 1
fi
grep -q "Use ModuleVelocityDistribution" "$mcc_file"
grep -q "PAPER_CASE_LEFT_BOUNDARY_RUNTIME" "$mcc_file"
grep -q "CF%NRun               =   $nt" "$app_dir/MCC_jw/input/controlflow.txt"
grep -q "PIC random seed" "$app_dir/code/PIC/Main_IFE_Test_2.f90"

cat <<EOF_DONE
Prepared standard paper baseline case:
  app_dir        = $app_dir
  distribution   = ${distribution}
  parameter      = ${distribution_shape}
  left_boundary  = $left_boundary
  species        = electron + Maxwellian ion (mi/me=${mass_ratio})
  ppc            = $ppg
  seed           = $random_seed
  dx, dt         = $dx, $dt
  nt             = $nt
  run mode       = ${run_mode}
  auto archive   = ${archive_results}
  Slurm CPUs     = 1
  Slurm memory   = ${requested_memory_mib} MiB
  Slurm walltime = ${requested_walltime}
  archive root   = $archive_root
  config         = $app_dir/case_config.txt
  slurm          = $app_dir/run_${case_name}.slurm

Expected immediate output:
  run.log
  run_metadata.txt
  OUTPUT/global_diagnostics.csv
EOF_DONE

if [[ "$archive_results" == "yes" ]]; then
  cat <<EOF_PRODUCTION

Expected paper snapshots and automatic archive:
  OUTPUT/Field/Average_x_${target_t30_file_step}.dat
  OUTPUT/Field/Average_x_${target_t50_file_step}.dat
  OUTPUT/Velocity/velocity_IJ_3${target_t30_file_step}.dat
  OUTPUT/Velocity/velocity_IJ_3${target_t50_file_step}.dat
  OUTPUT/Energy/energy_IJ2_${target_t30_file_step}.dat
  OUTPUT/Energy/energy_IJ2_${target_t50_file_step}.dat
  ${archive_root}/${case_name}
EOF_PRODUCTION
else
  cat <<EOF_SMOKE

This is a short smoke test. It verifies initialization, boundary sampling,
solver startup, and clean termination only; it is not a paper result.
EOF_SMOKE
fi
