using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;

namespace ParaxialWaveOptics {
	//レンズ
	public class Lens {
		public double distance;
		public double curvature;
		public double diameter;
		public double thickness;
		public double refractiveIndex;
		public bool isCurvatureRight; // 曲面が左右どちらか
		public double radius { get { return diameter / 2; } }
		public double focusLength { get { return curvature / (refractiveIndex - LensAnalyze.airRefractiveIndex); } }
	}
	public class Origin {
		public double position;
		public double height;
		public double divergence;
		public double wavelength;
	}
	//光源を(0,0)とした位置
	public class Position {
		public Position(double _x, double _y) {
			x = _x;
			y = _y;
		}
		public double x;
		public double y;
		public static Position operator +(Position p1, Position p2) {
			return new Position(p1.x + p2.x, p1.y + p2.y);
		}
		public static Position operator -(Position p1, Position p2) {
			return new Position(p1.x - p2.x, p1.y - p2.y);
		}
		public static Position operator *(Position p1, double coef) {
			return new Position(p1.x * coef, p1.y * coef);
		}
		public double length() {
			return Math.Sqrt(x * x + y * y);
		}
	}
	//入力
	public class ParameterSet {
		public Origin origin;
		public List<Lens> Lenss = new List<Lens>();
		public List<double> theta = new List<double>(); // 光源から追跡するレイの角度
		public List<double> time = new List<double>(); // 波面を追跡する時刻
	}

	//出力
	public class ResultSet {
		public List<List<Position>> wavefront = new List<List<Position>>(); // y方向の波面
		public List<List<Position>> raytrace = new List<List<Position>>(); // x方向のレイの軌跡
		public double lastAbberationWidth=0;
		public double lastAbberationHeight=0;
	}


	class LensAnalyze {
		//定数
		public const double vc = 3.0 * 1e8;
		public const double pi = Math.PI;
		public const double inch = 2.54 * 1e-2;
		public const double quartzRefractiveIndex = 1.53;
		public const double airRefractiveIndex = 1.0;
		public const double inputLensDiameter = inch * 2;

		private Position solveCrossPoint(double la,double lb,double lc,double ld,double centerx,double curvature, bool isRight) {
			if (isRight) {
				return new Position(
				((-la * ld) - lb * Math.Sqrt((la * la + lb * lb) * Math.Pow(curvature, 2) - ld * ld)) / (la * la + lb * lb) + centerx,
				((-lb * ld) + la * Math.Sqrt((la * la + lb * lb) * Math.Pow(curvature, 2) - ld * ld)) / (la * la + lb * lb)
				);
			} else {
				return new Position(
					((-la * ld) + lb * Math.Sqrt((la * la + lb * lb) * Math.Pow(curvature, 2) - ld * ld)) / (la * la + lb * lb) + centerx,
					((-lb * ld) - la * Math.Sqrt((la * la + lb * lb) * Math.Pow(curvature, 2) - ld * ld)) / (la * la + lb * lb)
					);
			}
		}
		public static int compareLens(Lens l1,Lens l2){
			if (l1.distance > l2.distance) return 1;
			else if (l1.distance < l2.distance) return -1;
			else if (l1.isCurvatureRight && !l2.isCurvatureRight) return 1;
			else if (!l1.isCurvatureRight && l2.isCurvatureRight) return -1;
			else return 0;
		}
		public ResultSet solveSimplest(ParameterSet ps) {
			ResultSet rs = new ResultSet();
			ps.Lenss.Sort(compareLens);
			for (int i = 0; i < ps.theta.Count; i++) {

				List<Position> wave = new List<Position>();
				List<Position> ray = new List<Position>();

				// first point
				int LensCount = 0;
				double currentTheta = ps.theta[i];
				double lastTime = 0;
				Position lastPosition = new Position(ps.origin.position, ps.origin.height);
				int timeIndex = 0;
				double opticalLength = 0;
				Position currentPosition = new Position(ps.origin.position, ps.origin.height);
				double currentTime = 0;
				ray.Add(currentPosition);

				while (timeIndex < ps.time.Count) {

					double distance = ps.Lenss[LensCount].distance-currentPosition.x;
					currentPosition = new Position(
						ps.Lenss[LensCount].distance,
						currentPosition.y + distance * currentTheta);
					currentTheta -= currentPosition.y / ps.Lenss[LensCount].focusLength;

					opticalLength = (currentPosition - lastPosition).length();
					currentTime += opticalLength / vc;

					while (timeIndex < ps.time.Count) {
						if (currentTime < ps.time[timeIndex]) break;
						wave.Add(lastPosition + (currentPosition - lastPosition) * ((ps.time[timeIndex] - lastTime) / (currentTime - lastTime)));
						timeIndex++;
					}

					lastTime = currentTime;
					lastPosition = currentPosition;
					ray.Add(currentPosition);
					LensCount++;

					if (LensCount >= ps.Lenss.Count) {
						// rest points
						while (timeIndex < ps.time.Count) {
							Position newPosition = new Position(
								currentPosition.x + vc * (ps.time[timeIndex] - currentTime) * Math.Cos(currentTheta),
								currentPosition.y + vc * (ps.time[timeIndex] - currentTime) * Math.Sin(currentTheta)
								);
							wave.Add(newPosition);
							if (timeIndex + 1 == ps.time.Count) ray.Add(newPosition);
							timeIndex++;
						}
					}
				}
				rs.wavefront.Add(wave);
				rs.raytrace.Add(ray);
			}


			return rs;
		}
		public ResultSet solveAbberation(ParameterSet ps) {

			ResultSet rs = new ResultSet();

			ps.Lenss.Sort(compareLens);

			for (int i = 0; i < ps.theta.Count; i++) {

				List<Position> wave = new List<Position>();
				List<Position> ray = new List<Position>();

				// first point
				bool isInside = false; // 屈折媒質の中にいるかどうか
				int LensCount = 0;
				double currentTheta = ps.theta[i];
				double lastTime = 0;
				Position lastPosition = new Position(ps.origin.position, ps.origin.height);
				int timeIndex = 0;
				double opticalLength = 0;
				Position currentPosition = new Position(ps.origin.position, ps.origin.height);
				double currentTime = 0;
				ray.Add(currentPosition);

				while (timeIndex < ps.time.Count) {

					bool isNextIntersectionCurvature = !(isInside ^ ps.Lenss[LensCount].isCurvatureRight);
					double relativeRefractiveIndex = ps.Lenss[LensCount].refractiveIndex / airRefractiveIndex;

					// 次の交点が円と直線の交点
					if (isNextIntersectionCurvature) {
						double centerx;
						double la, lb, lc, ld;
						bool isCurvaturePositive = ps.Lenss[LensCount].curvature > 0;
						// centerx 円の中心x座標
						// la x + lb y + lc = 0
						// ld = la centerx + lc

						centerx = ps.Lenss[LensCount].distance
							+ (isInside ? 1 : -1)*(ps.Lenss[LensCount].thickness
							- (isCurvaturePositive?1:-1) * Math.Sqrt(Math.Pow(ps.Lenss[LensCount].curvature, 2) - Math.Pow(ps.Lenss[LensCount].radius, 2)));

						la = Math.Tan(currentTheta);
						lb = -1.0;
						lc = lastPosition.y - Math.Tan(currentTheta) * lastPosition.x;
						ld = la * centerx + lc;
						currentPosition = solveCrossPoint(la, lb, lc, ld, centerx, ps.Lenss[LensCount].curvature,!(isInside ^ isCurvaturePositive));
		
						double lt = Math.Asin(currentPosition.y / Math.Abs(ps.Lenss[LensCount].curvature));
						if (isInside && isCurvaturePositive) 
							currentTheta = lt - Math.Asin(Math.Sin(lt - currentTheta) * relativeRefractiveIndex);
						if (!isInside && !isCurvaturePositive)
							currentTheta = lt - Math.Asin(Math.Sin(lt - currentTheta) / relativeRefractiveIndex);
						if (isInside && !isCurvaturePositive)
							currentTheta = Math.Asin(Math.Sin(currentTheta + lt) * relativeRefractiveIndex) - lt;
						if (!isInside && isCurvaturePositive)
							currentTheta = Math.Asin(Math.Sin(currentTheta + lt) / relativeRefractiveIndex) - lt;


//						if (isInside && !isCurvaturePositive) currentTheta = Math.Asin(Math.Sin(currentTheta + lt) / ps.Lenss[LensCount].refraxtionIndex) - lt;
//						else currentTheta = lt - Math.Asin(Math.Sin(lt - currentTheta) * ps.Lenss[LensCount].refraxtionIndex);

					}
					// 次の交点がy軸と直線の交点
					else {
						currentPosition = new Position(
							ps.Lenss[LensCount].distance,
							lastPosition.y + Math.Tan(currentTheta) * (ps.Lenss[LensCount].distance - lastPosition.x)
							);
						currentTheta = Math.Asin(Math.Sin(currentTheta) * (isInside ? relativeRefractiveIndex : 1.0/relativeRefractiveIndex));
					}

					if (Math.Abs(currentPosition.y) > ps.Lenss[LensCount].radius) currentPosition.y = double.NaN;

					opticalLength = (currentPosition - lastPosition).length();
					currentTime += opticalLength / vc * (isInside ? ps.Lenss[LensCount].refractiveIndex : airRefractiveIndex);

					while (timeIndex < ps.time.Count) {
						if (currentTime < ps.time[timeIndex]) break;
						wave.Add(lastPosition + (currentPosition - lastPosition) * ((ps.time[timeIndex] - lastTime) / (currentTime - lastTime)));
						timeIndex++;
					}

					lastTime = currentTime;
					lastPosition = currentPosition;
					ray.Add(currentPosition);

					if (isInside) LensCount++;
					isInside = !isInside;

					if (LensCount >= ps.Lenss.Count) {
						// rest points
						while (timeIndex < ps.time.Count) {
							Position newPosition = new Position(
								currentPosition.x + vc * (ps.time[timeIndex] - currentTime) * Math.Cos(currentTheta),
								currentPosition.y + vc * (ps.time[timeIndex] - currentTime) * Math.Sin(currentTheta)
								);
							wave.Add(newPosition);
							if (timeIndex + 1 == ps.time.Count) ray.Add(newPosition);
							timeIndex++;
						}
					}
				}
				rs.wavefront.Add(wave);
				rs.raytrace.Add(ray);
			}

			double maxZeroCross=double.MinValue;
			double minZeroCross=double.MaxValue;
			double absMinAngle = double.MaxValue;
			double absMinAngleZeroCross = 0;

			try {
				for (int i = 0; i < rs.raytrace.Count; i++) {
					if (ps.theta[i] == 0) continue;
					if (rs.raytrace[i].Count >= 2) {
						Position last = rs.raytrace[i][rs.raytrace[i].Count - 1];
						Position previous = rs.raytrace[i][rs.raytrace[i].Count - 2];
						if (last.y * previous.y <= 0) {
							double ratio = Math.Abs((previous.y - last.y) / previous.y);
							double zeroCross = previous.x + (last.x - previous.x) * ratio;
							maxZeroCross = Math.Max(zeroCross, maxZeroCross);
							minZeroCross = Math.Min(zeroCross, minZeroCross);

							if (absMinAngle > Math.Abs(ps.theta[i])) {
								absMinAngle = Math.Abs(ps.theta[i]);
								absMinAngleZeroCross = zeroCross;
							}

						} else throw new Exception();
					} else throw new Exception();
				}
				rs.lastAbberationWidth = maxZeroCross - minZeroCross;

				rs.lastAbberationHeight = double.MinValue;
				for (int i = 0; i < rs.raytrace.Count; i++) {
					if (ps.theta[i] == 0) continue;
					if (rs.raytrace[i].Count >= 2) {
						Position last = rs.raytrace[i][rs.raytrace[i].Count - 1];
						Position previous = rs.raytrace[i][rs.raytrace[i].Count - 2];
						if (last.y * previous.y <= 0) {
							double ratio = Math.Abs((last.x - previous.x) / (absMinAngleZeroCross-previous.x));
							double abberationHeight = previous.y + (last.y - previous.y) * ratio;

							rs.lastAbberationHeight = Math.Max(abberationHeight,rs.lastAbberationHeight);

						} else throw new Exception();
					} else throw new Exception();
				}
			} catch (Exception) {
				rs.lastAbberationWidth = double.MaxValue;
				rs.lastAbberationHeight = double.MaxValue;
			}
			return rs;
		}
	}
}
