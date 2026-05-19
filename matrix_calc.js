const EPS = 1e-8; // ค่าคลาดเคลื่อนเล็ก ๆ สำหรับกันหารด้วยศูนย์และเทียบค่าทศนิยม

// clone(A): คืนเมทริกซ์ใหม่ที่คัดลอกค่าจาก A เพื่อไม่แก้ข้อมูลต้นฉบับ
function clone(A){ return A.map(r=>[...r]); }

// eye(n): สร้างเมทริกซ์เอกลักษณ์ขนาด n x n
function eye(n){ return Array.from({length:n},(_,i)=>Array.from({length:n},(_,j)=>i===j?1:0)); }

// mul(A,B): คูณเมทริกซ์ A(n x k) กับ B(k x m) แล้วคืนผลลัพธ์ C(n x m)
function mul(A,B){
  const n=A.length,m=B[0].length,k=B.length;
  const C=Array.from({length:n},()=>Array(m).fill(0));
  for(let i=0;i<n;i++) for(let j=0;j<m;j++) for(let t=0;t<k;t++) C[i][j]+=A[i][t]*B[t][j];
  return C;
}

// transpose(A): สลับแถวเป็นคอลัมน์
function transpose(A){ return A[0].map((_,j)=>A.map(r=>r[j])); }

// norm(v): หาความยาวเวกเตอร์แบบ Euclidean = sqrt(sum(v_i^2))
function norm(v){ return Math.sqrt(v.reduce((s,x)=>s+x*x,0)); }

// qrDecompose(A): แตกเมทริกซ์ A เป็น Q และ R ด้วย Gram-Schmidt
function qrDecompose(A){
  const n=A.length, Q=Array.from({length:n},()=>Array(n).fill(0)), R=Array.from({length:n},()=>Array(n).fill(0));
  const V=transpose(A), U=[]; // V[i] คือคอลัมน์ที่ i ของ A
  for(let i=0;i<n;i++){
    let u=[...V[i]]; // เริ่มจากคอลัมน์เดิม
    for(let j=0;j<i;j++){
      const dot=V[i].reduce((s,x,k)=>s+x*U[j][k],0); // โปรเจกชันของ V[i] บน U[j]
      R[j][i]=dot;
      u=u.map((x,k)=>x-dot*U[j][k]); // ลบองค์ประกอบที่ซ้ำกับแกนเก่าออก
    }
    const un=norm(u)||EPS; // กันกรณี norm ใกล้ศูนย์
    R[i][i]=un;
    U.push(u.map(x=>x/un)); // normalize เป็นเวกเตอร์หน่วย
  }
  for(let i=0;i<n;i++) for(let j=0;j<n;j++) Q[j][i]=U[i][j]; // แปลงกลับเป็นรูปเมทริกซ์ Q
  return {Q,R};
}

// qrEigenvalues(A): ประมาณ eigenvalues โดยวน Ak+1 = RQ หลายรอบจนค่าแนวทแยงนิ่ง
function qrEigenvalues(A, it=220){
  let Ak=clone(A);
  for(let i=0;i<it;i++){
    const {Q,R}=qrDecompose(Ak);
    Ak=mul(R,Q);
  }
  return Ak.map((r,i)=>r[i]); // ค่าแนวทแยงของ Ak คือ eigenvalues โดยประมาณ
}

// rref(M): แปลงเมทริกซ์เป็น Reduced Row Echelon Form และคืนตำแหน่ง pivot
function rref(M){
  const A=clone(M); const rows=A.length, cols=A[0].length;
  let r=0; const pivots=[];
  for(let c=0;c<cols && r<rows;c++){
    let p=r;
    for(let i=r+1;i<rows;i++) if(Math.abs(A[i][c])>Math.abs(A[p][c])) p=i; // partial pivoting
    if(Math.abs(A[p][c])<EPS) continue; // คอลัมน์นี้ไม่มี pivot ใช้ไม่ได้
    [A[r],A[p]]=[A[p],A[r]]; // สลับแถวขึ้นมาเป็น pivot row
    const div=A[r][c]; for(let j=c;j<cols;j++) A[r][j]/=div; // ทำ pivot ให้เป็น 1
    for(let i=0;i<rows;i++) if(i!==r){
      const f=A[i][c]; if(Math.abs(f)<EPS) continue;
      for(let j=c;j<cols;j++) A[i][j]-=f*A[r][j]; // กำจัดค่าคอลัมน์ c ในแถวอื่นให้เป็น 0
    }
    pivots.push(c); r++;
  }
  return {A,pivots};
}

// eigenvectorsForLambda(A,l): หา basis ของ null space ของ (A - λI) => ได้ eigenvectors ของ λ
function eigenvectorsForLambda(A,l){
  const n=A.length;
  const M=A.map((row,i)=>row.map((v,j)=>v-(i===j?l:0))); // สร้าง A-λI
  const {A:R,pivots}=rref(M);
  const free=[]; for(let c=0;c<n;c++) if(!pivots.includes(c)) free.push(c); // ตัวแปรอิสระ
  const vecs=[];
  for(const f of free){
    const x=Array(n).fill(0); x[f]=1; // ตั้งตัวแปรอิสระทีละตัวเป็น 1
    for(let i=pivots.length-1;i>=0;i--){
      const pc=pivots[i]; let s=0;
      for(let j=pc+1;j<n;j++) s+=R[i][j]*x[j];
      x[pc]=-s; // back-substitution หา pivot variable
    }
    vecs.push(x);
  }
  return vecs;
}

// inverse(A): หา A^-1 ด้วยการทำ RREF บน [A | I]
function inverse(A){
  const n=A.length; const M=A.map((row,i)=>[...row,...eye(n)[i]]);
  const {A:R,pivots}=rref(M);
  if(pivots.length<n) return null; // rank ไม่เต็ม => อินเวิร์สไม่ได้
  return R.map(row=>row.slice(n)); // ด้านขวาหลัง reduce คือ A^-1
}

// uniqueApprox(vals): รวมค่า eigenvalue ที่เกือบเท่ากันให้เหลือค่าเดียว
function uniqueApprox(vals){
  const out=[];
  for(const v of vals) if(!out.some(x=>Math.abs(x-v)<1e-5)) out.push(v);
  return out;
}

// fmtMat(A): จัดรูปแบบเมทริกซ์ให้อ่านง่ายในรายงาน
function fmtMat(A){ return A.map(r=>`[ ${r.map(v=>v.toFixed(6)).join(' , ')} ]`).join('\n'); }

// buildReport(A): ฟังก์ชันหลัก สร้างรายงานผล diagonalization แบบข้อความ
export function buildReport(A){
  const n=A.length;
  if(![2,3].includes(n) || A.some(r=>r.length!==n)) throw new Error('Matrix must be 2x2 or 3x3');
  const vals = qrEigenvalues(A);
  const uniq = uniqueApprox(vals);

  let report = `A =\n${fmtMat(A)}\n\n1) Eigenvalues (approx)\n`;
  report += vals.map((v,i)=>`λ${i+1} = ${v.toFixed(6)}`).join('\n') + '\n\n2) Eigenvectors\n';

  const eigenPairs=[];
  for(const l of uniq){
    const vecs=eigenvectorsForLambda(A,l);
    report += `For λ = ${l.toFixed(6)}\n`;
    if(!vecs.length) report += '  - Could not compute a real eigenvector numerically.\n';
    vecs.forEach((v,idx)=>{ report += `  v${idx+1} = [${v.map(x=>x.toFixed(6)).join(', ')}]\n`; eigenPairs.push({l,v}); });
  }

  report += '\n3) Diagonalizable?\n';
  if(eigenPairs.length < n){
    report += 'No\n\n4) Cannot construct P and P^-1\n5) Reason: number of independent eigenvectors is less than matrix size.\n';
    return report;
  }

  const P = Array.from({length:n},(_,r)=>eigenPairs.slice(0,n).map(p=>p.v[r])); // ใช้ eigenvectors เป็นคอลัมน์ของ P
  const PInv = inverse(P);
  if(!PInv){
    report += 'No\n\n4) Cannot construct P and P^-1\n5) Reason: eigenvector matrix P is singular.\n';
    return report;
  }

  const D = mul(PInv,mul(A,P)); // ถ้า diagonalize ได้ D จะใกล้เมทริกซ์ทแยงมุม
  report += 'Yes\n\n4) P, P^-1, D = P^-1 A P\n';
  report += `P =\n${fmtMat(P)}\n\nP^-1 =\n${fmtMat(PInv)}\n\nD =\n${fmtMat(D)}\n\n5) If not diagonalizable, reason is shown above.\n`;
  return report;
}
