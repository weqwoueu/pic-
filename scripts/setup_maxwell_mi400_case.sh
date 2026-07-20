#!/usr/bin/env bash
set -euo pipefail

# Prepare the paper baseline plasma-expansion case:
#   distribution  : Maxwellian
#   species       : electron + one ion species
#   mass ratio    : mi/me = 400, Z = 1
#   domain        : [0, 1024 lambda_D0] x [0, 4 lambda_D0]
#   initial slab  : [0, 128 lambda_D0] x [0, 4 lambda_D0]
#   mesh          : dx = dy = 1.0 or 0.5 lambda_D0
#   time step     : dt = 0.05 or 0.025 / omega_pe
#
# Usage:
#   bash scripts/setup_maxwell_mi400_case.sh [ppc] [nt] [boundary] [seed] [dt] [dx]
#
# Examples:
#   bash scripts/setup_maxwell_mi400_case.sh 1000 20000 thermal 101 0.05 1.0
#   bash scripts/setup_maxwell_mi400_case.sh 80000 40000 thermal 101 0.025 1.0
#
# The paper-level ppc=80000 case is large. Use 1000 first as a standard
# configuration smoke test, then increase to 80000 for production.

ppg="${1:-1000}"
nt="${2:-20000}"
left_boundary="${3:-thermal}"
random_seed="${4:-101}"
dt="${5:-0.05}"
dx="${6:-1.0}"

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

if ! [[ "$ppg" =~ ^[1-9][0-9]*$ && "$nt" =~ ^[1-9][0-9]*$ && "$random_seed" =~ ^[0-9]+$ ]]; then
  echo "ppc and nt must be positive integers; seed must be a non-negative integer" >&2
  exit 2
fi

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
app_dir="$repo_root/PIC-IFE_GEC"
mcc_file="$app_dir/MCC_jw/code/Interface_IFE/MCCInterface.f90"

species_count=2
initial_particles=$((ppg * slab_cells_x * slab_cells_y * species_count))
particle_capacity=$(((initial_particles * 12 + 9) / 10))
if (( particle_capacity < 10000000 )); then
  particle_capacity=10000000
fi

target_t30_step="$(awk -v dt="$dt" 'BEGIN { printf "%.0f", 600.0 / dt }')"
target_t50_step="$(awk -v dt="$dt" 'BEGIN { printf "%.0f", 1000.0 / dt }')"
printf -v target_t30_file_step '%06d' "$target_t30_step"
printf -v target_t50_file_step '%06d' "$target_t50_step"

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

# Maxwellian initialization.
perl -0pi -e 's/Integer\(4\) :: MaxKappa=[0-9]+/Integer(4) :: MaxKappa=1/' "$mcc_file"

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

if [[ "$left_boundary" == "thermal" ]]; then
  left_block='        If (PO%X <= dxmin) Then
            ! PAPER_CASE_LEFT_BOUNDARY_BEGIN thermal-reservoir Maxwellian
            ek_before = PO%Energy(PB%Mass, PB%VFactor)
            Iposflag = 1
            TimeRemain = (PO%X - dxmin)/PO%Vx
            TimeMove = TimeRemain
            PO%Y = PO%Y - PO%Vy * TimeRemain
            PO%Z = PO%Z - PO%Vz * TimeRemain
            PO%X = dxmin + 10E-6
            call DRandom(ranf1)
            call DRandom(ranf2)
            call DRandom(ranf3)
            PO%Vx = SQRT((-DLOG(Max(ranf1, 1.0d-300)))/beta)
            velocity_tangential = SQRT((-DLOG(Max(ranf2, 1.0d-300)))/beta)
            theta_v = 2*pii*ranf3
            PO%Vy = velocity_tangential*COS(theta_v)
            PO%Vz = velocity_tangential*SIN(theta_v)
            VFactor = 1.0 / PB%VFactor
            call PO%VelRes(VFactor)
            ek_after = PO%Energy(PB%Mass, PB%VFactor)
            If (PB%UnequalWeightFlag) Then
                diag_weight = PO%WQ
            Else
                diag_weight = PB%Weight
            End If
            Call AddDiagThermalExchange(diag_weight, ek_before, ek_after)
            ! PAPER_CASE_LEFT_BOUNDARY_END'
else
  left_block='        If (PO%X <= dxmin) Then
            ! PAPER_CASE_LEFT_BOUNDARY_BEGIN specular
            Iposflag = 1
            TimeRemain = (PO%X - dxmin)/PO%Vx
            TimeMove = TimeRemain
            PO%Y = PO%Y - PO%Vy * TimeRemain
            PO%Z = PO%Z - PO%Vz * TimeRemain
            PO%Vx = -PO%Vx
            PO%X = dxmin + 10E-6
            ! PAPER_CASE_LEFT_BOUNDARY_END'
fi
export PAPER_LEFT_BOUNDARY_BLOCK="$left_block"
replace_block "$mcc_file" 'BEGIN { $r = $ENV{"PAPER_LEFT_BOUNDARY_BLOCK"} . "\n" } s/        If \(PO%X <= dxmin\) Then.*?        Elseif \(PO%X > dxmax\) Then/$r        Elseif (PO%X > dxmax) Then/s'
unset PAPER_LEFT_BOUNDARY_BLOCK

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

ion_mass="3.6438D-28"
cat > "$app_dir/INPUT/pic.inp" <<EOF_PIC
! Maxwellian mi/me=400 two-species plasma expansion
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

case_name="maxwellian_dx${dx_tag}_dt${dt_tag}_ppc${ppg}_seed${random_seed}_${left_boundary}"

cat > "$app_dir/case_config.txt" <<EOF_CONFIG
case_name = ${case_name}
distribution = maxwellian
left_boundary = ${left_boundary}
species = electron, ion
ion_mass_over_electron_mass = 400
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
particles_per_cell_per_species = ${ppg}
initial_particles_total = ${initial_particles}
particle_capacity = ${particle_capacity}
target_omega_pi_t_30_step = ${target_t30_step}
target_omega_pi_t_50_step = ${target_t50_step}
diagnostic_stride_global = 1000
field_velocity_output_stride = 1000
random_seed = ${random_seed}
EOF_CONFIG

cat > "$app_dir/run_${case_name}.slurm" <<EOF_SLURM
#!/bin/bash
#SBATCH -J mx_${dx_tag}_${dt_tag}_${ppg}
#SBATCH -p comp
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH -t 36:00:00
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.err

module purge
module load intel/2022.1
module load cmake/3.23.5

cd /data/home/dg001947/pic-/PIC-IFE_GEC
ulimit -s unlimited
export PIC_RANDOM_SEED=${random_seed}

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
git_revision="\$(git -C /data/home/dg001947/pic- rev-parse HEAD 2>/dev/null || echo unknown)"

cat > run_metadata.txt <<EOF_METADATA
case_name = ${case_name}
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

cd /data/home/dg001947/pic-
bash scripts/archive_verification_case.sh
EOF_SLURM

pic_ispe="$(awk 'NR==11 {print $1}' "$app_dir/INPUT/pic.inp" | tr -d ',')"
pic_nt="$(awk 'NR==12 {print $1}' "$app_dir/INPUT/pic.inp" | tr -d ',')"
if [[ "$pic_ispe" != "2" || "$pic_nt" != "$nt" ]]; then
  echo "case validation failed: ispe_tot=$pic_ispe nt=$pic_nt" >&2
  exit 1
fi
grep -q "Integer(4) :: MaxKappa=1" "$mcc_file"
grep -q "PAPER_CASE_LEFT_BOUNDARY_BEGIN" "$mcc_file"
grep -q "CF%NRun               =   $nt" "$app_dir/MCC_jw/input/controlflow.txt"
grep -q "PIC random seed" "$app_dir/code/PIC/Main_IFE_Test_2.f90"

cat <<EOF_DONE
Prepared standard paper baseline case:
  app_dir        = $app_dir
  distribution   = Maxwellian
  left_boundary  = $left_boundary
  species        = electron + ion (mi/me=400)
  ppc            = $ppg
  seed           = $random_seed
  dx, dt         = $dx, $dt
  nt             = $nt
  config         = $app_dir/case_config.txt
  slurm          = $app_dir/run_${case_name}.slurm

Expected key outputs:
  OUTPUT/global_diagnostics.csv
  OUTPUT/Field/Average_x_${target_t30_file_step}.dat
  OUTPUT/Field/Average_x_${target_t50_file_step}.dat
  OUTPUT/Velocity/velocity_IJ_3${target_t30_file_step}.dat
  OUTPUT/Velocity/velocity_IJ_3${target_t50_file_step}.dat
EOF_DONE
