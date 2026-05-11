#!/usr/bin/env bash
set -euo pipefail

# Prepare the collisionless plasma-expansion case from the paper:
# Maxwellian initial velocity, specular/full reflection on the left x boundary,
# outflow on the right x boundary, y-periodic particles, and mi/me = 400.
#
# Usage:
#   bash scripts/setup_maxwell_mi400_case.sh [particles_per_cell] [nt]
#
# For a quick proof run, use 200-1000 particles/cell. The paper value is 80000.

ppg="${1:-1000}"
nt="${2:-20000}"

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
app_dir="$repo_root/PIC-IFE_GEC"

backup_once() {
  local f="$1"
  if [[ -f "$f" && ! -f "$f.maxwell_mi400.bak" ]]; then
    cp "$f" "$f.maxwell_mi400.bak"
  fi
}

backup_once "$app_dir/code/PIC/Main_IFE_Test_2.f90"
backup_once "$app_dir/MCC_jw/code/Interface_IFE/MCCInterface.f90"
backup_once "$app_dir/code/In-Output/Output_To_Tecplot_IJK_2D.f90"
backup_once "$app_dir/OUTPUT_velocity.f90"
backup_once "$app_dir/Output_Energy.f90"
backup_once "$app_dir/INPUT/mesh.inp"
backup_once "$app_dir/INPUT/object.inp"
backup_once "$app_dir/INPUT/pic.inp"
backup_once "$app_dir/MCC_jw/input/gas.txt"
backup_once "$app_dir/MCC_jw/input/controlflow.txt"

# Initial slab [0,128] x [0,4], embedded in the full [0,1024] x [0,4] domain.
perl -0pi -e 's/dxmaxmax = 200\.0/dxmaxmax = 128.0/' "$app_dir/code/PIC/Main_IFE_Test_2.f90"
perl -0pi -e 's/dymaxmax = 4\.0/dymaxmax = ymax/' "$app_dir/code/PIC/Main_IFE_Test_2.f90"
perl -0pi -e 's/    call OUTPUT_velocity\(it\)\n    call Output_Energy\(it\)/    If (Mod(it,1000).eq.0 .Or. it == 1 .Or. it == nt) Then\n        call OUTPUT_velocity(it)\n        call Output_Energy(it)\n    Endif/' "$app_dir/code/PIC/Main_IFE_Test_2.f90"
perl -0pi -e 's/    If \(Mod\(it,1000\)\.eq\.0 \.Or\. it == 1\) Then/    If (Mod(it,1000).eq.0 .Or. it == 1 .Or. it == nt) Then/' "$app_dir/code/PIC/Main_IFE_Test_2.f90"

# Maxwellian initialization.
perl -0pi -e 's/Integer\(4\) :: MaxKappa=2/Integer(4) :: MaxKappa=1/' "$app_dir/MCC_jw/code/Interface_IFE/MCCInterface.f90"

# Make the MCC/JW particle bundle use the species data read from pic.inp.
perl -0pi -e 's/    Use Field_2D, Only:dens0/    Use Field_2D, Only: dens0/' "$app_dir/MCC_jw/code/Interface_IFE/MCCInterface.f90"
perl -0pi -e 's/\n    Use Particle_2D, Only: qs, xm, tmpj, ispe_tot//' "$app_dir/MCC_jw/code/Interface_IFE/MCCInterface.f90"
if ! grep -q 'PIC_EXPANSION_MI400' "$app_dir/MCC_jw/code/Interface_IFE/MCCInterface.f90"; then
  perl -0pi -e 's/(        Call GasInitPegasus\(ControlFlowGlobal\).*\n)/$1        ! PIC_EXPANSION_MI400: use pic.inp masses, charges, and temperatures.\n        Do i = 0, Min(ControlFlowGlobal%Ns, ispe_tot - 1)\n            SpecyGlobal(i)%Charge = qs(i+1) * q_ref\n            SpecyGlobal(i)%Mass = xm(i+1) * m_ref\n            If (Allocated(tmpj)) Then\n                SpecyGlobal(i)%InitTemperature = tmpj(i+1) * T_ref\n                SpecyGlobal(i)%Temperature = SpecyGlobal(i)%InitTemperature\n            End If\n            If (Allocated(dens0)) Then\n                SpecyGlobal(i)%InitDensity = dens0(i+1)\n                SpecyGlobal(i)%Density = SpecyGlobal(i)%InitDensity\n            End If\n        End Do\n/' "$app_dir/MCC_jw/code/Interface_IFE/MCCInterface.f90"
fi

# Let ParticlePerGrid in pic.inp control the initial slab particle count.
perl -0pi -e 's/        !======= for wsy paper case ==========\n        PB%NParNormal = 6000000\n        NPArMax = PB%NParNormal\n        If\(Allocated\(PB%PO\)\) Deallocate\(PB%PO\)\n        Allocate\(PB%PO\(NPArMax\)\)\n        !======================================\n/        ! Initial particle storage is allocated below from ParticlePerGrid.\n/' "$app_dir/MCC_jw/code/Interface_IFE/MCCInterface.f90"
perl -0pi -e 's/           PB%NPar = 5000 \* \(dxmaxmax-dxminmin\) \* \(dymaxmax-dyminmin\)\s*\n            PB%Weight = affp_bjw\(isp\)\s*\n            !PB%Weight = dens0\(isp\)\*n_ref \* RegionVolume \/ PB%NPar\s*\n            !print\*,dens0\(isp\)\*n_ref\s*\n            PB%NParNormal = 5000 \* \(dxmaxmax-dxminmin\) \* \(dymaxmax-dyminmin\)[^\n]*\n/            If (delta_global == 0) Then\n                RegionVolume = (dxmaxmax-dxminmin)*L_ref * (dymaxmax-dyminmin)*L_ref\n            Elseif (delta_global == 1) Then\n                RegionVolume = (dxmaxmax-dxminmin)*L_ref * PI*(dymaxmax**2-dyminmin**2)*L_ref**2\n            Endif\n            PB%NPar = INT(DBLE(ParticlePerGrid) * (dxmaxmax-dxminmin) * (dymaxmax-dyminmin))\n            PB%Weight = dens0(isp)*n_ref * RegionVolume \/ DBLE(PB%NPar)\n            PB%NParNormal = PB%NPar\n/s' "$app_dir/MCC_jw/code/Interface_IFE/MCCInterface.f90"
perl -0pi -e 's/NPArMax = Ceiling\(3\.0 \* PB%NParNormal\)/NPArMax = Ceiling(1.2D0 * PB%NParNormal)/' "$app_dir/MCC_jw/code/Interface_IFE/MCCInterface.f90"

# Left x boundary: specular/full reflection. Right x boundary already deletes escaped particles.
perl -0pi -e 's/            Iposflag = 1\n            TimeRemain = \(PO%X - dxmin\)\/PO%Vx\n            TimeMove = TimeRemain\n            PO%Y = PO%Y - PO%Vy \* TimeRemain\n            PO%X = dxmin \+ 10E-6\n            !call PO%VelKappaInit[^\n]*\n            call DRandom\(ranf1\)\n            CALL RANDOM_NUMBER\(ranf1\)\n            PO%Vx =SQRT \(\(1\.d0\*Kappa-1\.5\)\/beta\*\(\(ranf1\)\*\*\(-1\/\(Kappa-1\)\)-1\)\)\n            VFactor = 1\.0 \/ PB%VFactor\n            call PO%VelRes\(VFactor\)/            Iposflag = 1\n            TimeRemain = (PO%X - dxmin)\/PO%Vx\n            TimeMove = TimeRemain\n            PO%Y = PO%Y - PO%Vy * TimeRemain\n            PO%Z = PO%Z - PO%Vz * TimeRemain\n            PO%Vx = -PO%Vx\n            PO%X = dxmin + 10E-6/' "$app_dir/MCC_jw/code/Interface_IFE/MCCInterface.f90"

# y boundaries: periodic particle wrap, matching the paper setup.
perl -0pi -e 's/        Elseif \(PO%Y < dymin\) Then.*?        Elseif \(PO%Y > dymax\) Then.*?        Endif\n        \n        If \(N_objects/        Elseif (PO%Y < dymin) Then\n            If (delta == 0) Then\n                Do While (PO%Y < dymin)\n                    PO%Y = PO%Y + (dymax - dymin)\n                End Do\n            Elseif (delta == 1) Then\n                Print*,'"'"'wrong cross'"'"'\n                Print*,PO\n                pause\n            Endif\n        Elseif (PO%Y > dymax) Then\n            If (delta == 0) Then\n                Do While (PO%Y > dymax)\n                    PO%Y = PO%Y - (dymax - dymin)\n                End Do\n            Endif\n        Endif\n        \n        If (N_objects/s' "$app_dir/MCC_jw/code/Interface_IFE/MCCInterface.f90"

# Two-species output: avoid the historical third-species columns.
perl -0pi -e 's/WRITE\(1,\('\''\(A150\)'\''\)\) '\''VARIABLES = "x" "y" "R" "Phi" "Rho" "Rho_ele" "Rho_H" "Rho_C" "efx" "efy" "Ek_ele" "Ek_ion" "Ek_tot_ele" "Ek_tot_ion" "node_type1"  "node_type2"'\''/WRITE(1,('\''(A150)'\'')) '\''VARIABLES = "x" "y" "R" "Phi" "Rho" "Rho_ele" "Rho_ion" "efx" "efy" "Ek_ele" "Ek_ion" "Ek_tot_ele" "Ek_tot_ion" "node_type1"  "node_type2"'\''/' "$app_dir/code/In-Output/Output_To_Tecplot_IJK_2D.f90"
perl -0pi -e 's/          WRITE\(1,50\) HP\(1:2,i\), radius\(i\),Phi\(i,1\), Rho\(i,1\), Rho_s\(i,1,1\), Rho_s\(i,1,2\),Rho_s\(i,1,3\), efx\(i, 1\), efy\(i, 1\), Ek_s\(i,1,1\), Ek_s\(i,1,2\), Ek_tot\(i,1,1\), Ek_tot\(i,1,2\), node_type\(1, i\), node_type\(2, i\)/          WRITE(1,*) HP(1:2,i), radius(i), Phi(i,1), Rho(i,1), Rho_s(i,1,1), Rho_s(i,1,2), \&\n                     efx(i,1), efy(i,1), Ek_s(i,1,1), Ek_s(i,1,2), Ek_tot(i,1,1), Ek_tot(i,1,2), \&\n                     node_type(1,i), node_type(2,i)/' "$app_dir/code/In-Output/Output_To_Tecplot_IJK_2D.f90"

perl -0pi -e 's/        do num=1,N\n            !WRITE\(544, '\''\(A150\)'\''\) [^\n]*\n            WRITE\(544, 18\) driftVelocity\(1,num, 1\), driftVelocity\(2,num, 1\), driftVelocity\(3,num, 1\), thermalVelocity\(1,num, 1\), thermalVelocity\(2,num, 1\), thermalVelocity\(3,num, 1\),x_num\(1,num\),x_num\(2,num\),x_num\(3,num\),num\n        end do\n    \n        \n18      FORMAT[^\n]*\n/        WRITE(544, '\''(A)'\'') '\''# drift_e drift_i thermal_e thermal_i count_e count_i cell'\''\n        do num=1,N\n            WRITE(544, *) driftVelocity(1,num,1), driftVelocity(2,num,1), thermalVelocity(1,num,1), \&\n                          thermalVelocity(2,num,1), x_num(1,num), x_num(2,num), num\n        end do\n    \n        \n/s' "$app_dir/OUTPUT_velocity.f90"
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
1025, 5, 0
1.0, 1.0
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
! Maxwellian specular-reflection plasma expansion, mi/me = 400
0, 0
.false.
.false.
1.0D21
1.0
$ppg
1.0D0, 1.0D0, 1.0D0
10000000
.FALSE., .FALSE.
2
$nt, 0.05
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
/
EOF_CF

mkdir -p "$app_dir/OUTPUT/Field" "$app_dir/OUTPUT/Velocity" "$app_dir/OUTPUT/Energy" "$app_dir/OUTPUT/Average" "$app_dir/DUMP"

cat <<EOF_DONE
Prepared Maxwellian mi/me=400 case in:
  $app_dir

Particles per cell: $ppg
nt: $nt

For the paper-level particle count, rerun:
  bash scripts/setup_maxwell_mi400_case.sh 80000 20000
EOF_DONE
