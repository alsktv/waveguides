# Coupled Waveguide Electromagnetics Simulator 🌊

이중 도파로 전자기파 시뮬레이터는 두 개의 평행한 유전체 평판 도파로(Slab Waveguide) 사이의 **에바네센트 파 결합(Evanescent Wave Coupling)** 및 이에 따른 전력 전이 현상을 실시간으로 시각화하고 물리학적으로 분석해주는 인터랙티브 웹 애플리케이션입니다.

본 시뮬레이터는 **Cloudflare Pages**를 통해 호스팅 및 배포되었습니다.

## 🌟 주요 기능

- **실시간 전기장 강도(Electric Field Intensity) 2D 렌더링**:
  - 도파로의 구조(코어/클래딩 경계면)를 투영하고, 진행하는 전자기파의 전파 양상을 60fps 속도로 렌더링합니다.
  - 전위 위상에 따라 양의 전기장은 하늘색/청색, 음의 전기장은 핑크색/자주색으로 직관적인 컬러 맵 매핑을 제공합니다.
- **다양한 광 입사 모드(Launch Mode) 지원**:
  - **도파로 1 입사 (WG1)** 또는 **도파로 2 입사 (WG2)**: 에너지가 도파로 간에 교차하여 Sinusoidal 형태로 왕복 이동하는 현상을 관찰할 수 있습니다.
  - **대칭 모드 (Symmetric Supermode / Even)** 및 **반대칭 모드 (Antisymmetric Supermode / Odd)**: 결합이 고정되어 에너지의 수평적 이동 없이 정상파 프로파일로 진행하는 양상을 나타냅니다.
- **동적 파라미터 튜닝**:
  - 도파로 간격(Gap), 도파로 너비(Width), 입사 광 파장($\lambda$), 클래딩과의 굴절률 차이($\Delta n$) 등의 조절 슬라이더를 제공하며 조작 시 실시간으로 파동 방정식의 해석적 결과가 즉시 연동됩니다.
- **물리적 신뢰도**:
  - 단일 유전체 평판 도파로의 초월 함수 방정식(Transcendental Equation)인 $u\tan(u) = \sqrt{V^2 - u^2}$을 Bisection 수치 해석으로 매 순간 풀어내어 횡방향 가이딩 파라미터($k_x$, $\gamma$, $\beta$)를 결정하고, 해석적 결합 모드 이론(Coupled Mode Theory) 수식을 통해 물리적으로 완벽한 결합 계수($\kappa$)를 계산합니다.
- **실시간 종방향 파워 그래프**:
  - 전자기파가 도파로를 따라 진행함에 따른 에너지 전력 분포($P_1(z)$, $P_2(z)$)의 파동 그래프를 하단에 연동하여 도시해 줍니다.

## 📂 파일 구조

- `index.html`: 레이아웃 뼈대 마크업 및 LaTeX 공식을 표시하기 위한 MathJax 수식 패키지 로드.
- `style.css`: 그래디언트 및 글래스모피즘(Glassmorphism) 효과를 가미한 대시보드 다크 테마 스타일시트.
- `script.js`: 유전체 도파로 파동 방정식을 계산하는 물리 엔진 및 Canvas 픽셀 버퍼 애니메이션 렌더링 스크립트.

## 🚀 Cloudflare Pages 배포 방법

이 프로젝트는 static HTML 프로젝트이므로, Cloudflare Pages와 연동하여 1분 이내에 빌드 및 배포할 수 있습니다.

### 🌐 Cloudflare Dashboard 방식 (추천)
1. [Cloudflare 대시보드](https://dash.cloudflare.com/)에 로그인합니다.
2. 사이드바에서 **Workers & Pages** 메뉴로 이동합니다.
3. **Create** -> **Pages** -> **Connect to Git** 버튼을 차례로 클릭합니다.
4. 본인의 GitHub 계정을 연동한 후, `waveguides` 저장소(`https://github.com/alsktv/waveguides.git`)를 선택합니다.
5. **Set up builds and deployments** 단계에서 다음과 같이 설정합니다:
   - **Framework preset**: `None`
   - **Build command**: (비워둠)
   - **Build output directory**: `.` (루트 경로)
6. **Save and Deploy** 버튼을 누르면 배포 파이프라인이 즉시 시작되어 고유한 `*.pages.dev` 도메인 주소로 시뮬레이터가 무료로 영구 호스팅됩니다!
