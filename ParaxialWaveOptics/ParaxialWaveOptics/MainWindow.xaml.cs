using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Data;
using System.Windows.Documents;
using System.Windows.Input;
using System.Windows.Media;
using System.Windows.Media.Imaging;
using System.Windows.Navigation;
using System.Windows.Shapes;
using System.Collections.ObjectModel;
using System.ComponentModel;
using System.Data;
using System.Windows.Controls.Primitives;


namespace ParaxialWaveOptics {
	/// <summary>
	/// MainWindow.xaml の相互作用ロジック
	/// </summary>
	public partial class MainWindow : Window, INotifyPropertyChanged {

		public event PropertyChangedEventHandler PropertyChanged;
		// Create the OnPropertyChanged method to raise the event
		protected void OnPropertyChanged(string name) {
			PropertyChangedEventHandler handler = PropertyChanged;
			if (handler != null) {
				handler(this, new PropertyChangedEventArgs(name));
			}
		}

		public int rayCount { get; set; }
		public int waveCount { get; set; }
		public string minDistanceString { get; set; }
		private double minDistance { get { return double.Parse(minDistanceString) * 1e-2; } }
		public string maxDistanceString { get; set; }
		private double maxDistance { get { return double.Parse(maxDistanceString) * 1e-2; } }
		public string drawHeightString { get; set; }
		private double drawHeight { get { return double.Parse(drawHeightString) * LensAnalyze.inch; } }

		double maxAngle, minAngle;
		public string abberationWidthString_;
		public string abberationWidthString { get { return abberationWidthString_; } set { abberationWidthString_ = value; OnPropertyChanged("abberationWidthString"); } }
		public string abberationHeightString_;
		public string abberationHeightString { get { return abberationHeightString_; } set { abberationHeightString_ = value; OnPropertyChanged("abberationHeightString"); } }
		public class LensText: INotifyPropertyChanged {
			public event PropertyChangedEventHandler PropertyChanged;
			// Create the OnPropertyChanged method to raise the event
			protected void OnPropertyChanged(string name) {
				PropertyChangedEventHandler handler = PropertyChanged;
				if (handler != null) {
					handler(this, new PropertyChangedEventArgs(name));
				}
			}
			public bool isEnabled { get; set; }
			public string distance_;
			public string distance { get { return distance_; } set { distance_ = value; OnPropertyChanged("distance"); } }
			public string curvature { get; set; }
			public string diameter { get; set; }
			public string refractiveIndex { get; set; }
			public bool isCurvatureRight { get; set; } // 曲面が左右どちらか
			public string thickness { get; set; }
		}
		public ObservableCollection<LensText> LensData;
		public class OriginText {
			public bool isEnabled { get; set; }
			public string position { get;set;}
			public string height { get; set; }
			public string divergence { get; set; }
			public string wavelength { get; set; }
			public string color { get; set; }
		}
		public ObservableCollection<OriginText> originData;
		public MainWindow() {
			InitializeComponent();

			LensData = new ObservableCollection<LensText>(){ 
				new LensText() { 
					isEnabled = true, 
					distance = "3", 
					curvature = "2.0", 
					thickness="5.0",
					diameter = "1" ,
					refractiveIndex = LensAnalyze.quartzRefractiveIndex.ToString(), 
					isCurvatureRight = true
				},
				new LensText() {
					isEnabled=true, 
					distance = "20", 
					curvature = "2.5",
					thickness="5.0",
					diameter = "1", 
					refractiveIndex = LensAnalyze.quartzRefractiveIndex.ToString(), 
					isCurvatureRight = false 
				}
			};
			originData = new ObservableCollection<OriginText>(){
				new OriginText(){
					isEnabled=true,
					position="0",
					height="0",
					divergence="0.05",
					wavelength="780",
					color="00EE00"
				},
				new OriginText(){
					isEnabled=true,
					position="0",
					height="2",
					divergence="0.05",
					wavelength="780",
					color="0000EE"
				}
			};

			rayCount = 4;
			waveCount = 4;
			minDistanceString = "0";
			maxDistanceString = "30";
			drawHeightString = "1.5";
			Text_RayCount.DataContext = this;
			Text_WaveCount.DataContext = this;
			Text_MaxDistance.DataContext = this;
			Text_MinDistance.DataContext = this;
			Text_DrawHeight.DataContext = this;
			Label_AbberationHeight.DataContext = this;
			Label_AbberationWidth.DataContext = this;

			LensList.ItemsSource = LensData;
			OriginList.ItemsSource = originData;
		}

		private void Callback_WindowResized(object sender, SizeChangedEventArgs e) {
			update();
		}
		private void Callback_RepaintButton(object sender, RoutedEventArgs e) {
			update();
		}
		private void Callback_ParameterChanged(object sender, RoutedEventArgs e) {
			update();
		}
		private void Callback_DataChanged(object sender, DataGridCellEditEndingEventArgs e) {
			update();
		}
		private void Callback_Update(object sender, DataTransferEventArgs e) {
			update();
		}
		private void update() {
			try {

				solve();

//				abberationHeightString = (res[0].lastAbberationHeight * 1e3).ToString("G3");
//				abberationWidthString = (res[0].lastAbberationWidth * 1e3).ToString("G3");

			} catch (Exception) {
				List<UIElement> deleted = new List<UIElement>();
				foreach (UIElement uie in drawCanvas.Children) {
					if (uie is Line || uie is Ellipse) deleted.Add(uie);
				}
				foreach (UIElement uie in deleted) {
					drawCanvas.Children.Remove(uie);
				}
			}
		}

		private List<ResultSet> solve() {
			ParameterSet param;
			param = new ParameterSet();
			for (int i = 0; i < LensData.Count; i++) {
				if (LensData[i].isEnabled) {
					Lens Lens = new Lens();
					Lens.distance = double.Parse(LensData[i].distance) * 1e-2;
					Lens.curvature = double.Parse(LensData[i].curvature) * 1e-2;
					Lens.diameter = double.Parse(LensData[i].diameter) * LensAnalyze.inch;
					Lens.thickness = double.Parse(LensData[i].thickness) * 1e-3;
					Lens.refractiveIndex = double.Parse(LensData[i].refractiveIndex);
					Lens.isCurvatureRight = LensData[i].isCurvatureRight;
					Lens.diameter = Math.Min(Lens.curvature * 2, Lens.diameter);
					param.Lenss.Add(Lens);
				}
			}
			if (param.Lenss.Count <= 0) throw new Exception();
			double opticalMaxDistance = maxDistance;
			double opticalMinDistance = minDistance;
			for (int i = 0; i < param.Lenss.Count;i++){
				double cuv = param.Lenss[i].curvature;
				double rad = param.Lenss[i].radius;
				double thi = param.Lenss[i].thickness;
				double pos = param.Lenss[i].distance;
				double dif = Math.Abs(cuv) - Math.Sqrt(cuv*cuv-rad*rad);
				double rf = param.Lenss[i].refractiveIndex;
				double left, right;
				if (param.Lenss[i].isCurvatureRight && param.Lenss[i].curvature > 0) {
					left = pos; right = pos + thi + dif;
				}else if (!param.Lenss[i].isCurvatureRight && param.Lenss[i].curvature > 0) {
					left = pos - dif - thi; right = pos;
				}else if (param.Lenss[i].isCurvatureRight && !(param.Lenss[i].curvature > 0)) {
					left = pos; right = pos + thi - dif;
				}else {
					left = pos - thi + dif; right = pos;
				}
				if (left < minDistance && minDistance < right) opticalMinDistance += (minDistance - left) * (rf - LensAnalyze.airRefractiveIndex);
				if (right < minDistance) opticalMinDistance += (right - left) * (rf - LensAnalyze.airRefractiveIndex);
				if (left < maxDistance && maxDistance < right) opticalMaxDistance += (maxDistance - left) * (rf - LensAnalyze.airRefractiveIndex);
				if (right < maxDistance) opticalMaxDistance += (right - left) * (rf - LensAnalyze.airRefractiveIndex);
			}
			for (int j = 0; j < waveCount; j++) {
				param.time.Add(opticalMaxDistance / LensAnalyze.vc + (opticalMaxDistance - opticalMinDistance) / LensAnalyze.vc * j / (waveCount - 1));
			}

			clearPaint();
			paintCommon(param);
			List<ResultSet> ress = new List<ResultSet>();
			for (int orgc = 0; orgc < originData.Count; orgc++) {
				if (originData[orgc].isEnabled) {
					param.origin = new Origin();
					param.origin.position = double.Parse(originData[orgc].position) * 1e-2;
					param.origin.height = double.Parse(originData[orgc].height) * 1e-3;
					param.origin.wavelength = double.Parse(originData[orgc].wavelength) * 1e-9;

					double divergence = double.Parse(originData[orgc].divergence);
					maxAngle = divergence;
					minAngle = -divergence;
//					maxAngle = Math.Atan((1.0 + atomHeight / LensAnalyze.inch) * Math.Tan(inputAngle));
//					minAngle = -Math.Atan((1.0 - atomHeight / LensAnalyze.inch) * Math.Tan(inputAngle));
//					maxAngle = Math.Atan((param.Lenss[0].radius-atomHeight) / param.Lenss[0].distance);
//					minAngle = Math.Atan((-param.Lenss[0].radius-atomHeight) / param.Lenss[0].distance);
					param.theta.Clear();
					for (int j = 0; j < rayCount; j++) {
						param.theta.Add(minAngle + (maxAngle - minAngle) * j / (rayCount - 1));
					}

					if (rayCount == 0 || waveCount == 0) continue;
//					double inputAngle = Math.Atan(LensAnalyze.inputLensDiameter / 2 / inputFocalLength);

					ResultSet res;
					if (CheckBox_CalcMethod.IsChecked.Value) {
						res = (new LensAnalyze()).solveAbberation(param);
					} else {
						res = (new LensAnalyze()).solveSimplest(param);
					}
					int colori = Convert.ToInt32(originData[orgc].color,16);
					byte r = (byte)(colori % 0x100);
					byte g = (byte)((colori / 0x100) % 0x100);
					byte b = (byte)((colori / 0x100 / 0x100) % 0x100);
					Color color = new Color { R = r, G = g, B = b };
					paintResult(param,res,color);
					ress.Add(res);

				}
			}
			return ress;
		}

		private void clearPaint() {
			List<UIElement> deleted = new List<UIElement>();
			foreach (UIElement uie in drawCanvas.Children) {
				if (uie is Line || uie is Ellipse) deleted.Add(uie);
			}
			foreach (UIElement uie in deleted) {
				drawCanvas.Children.Remove(uie);
			}
		}
		private void paintCommon(ParameterSet param) {
			double width = drawCanvas.ActualWidth;
			double height = drawCanvas.ActualHeight;

			double xcoef = width / (maxDistance - minDistance);
			double ycoef = height / drawHeight;
			double xoffset = -minDistance * xcoef;

			Line gridLine = new Line() { Stroke = Brushes.LightGray, StrokeThickness = 2 };
			gridLine.X1 = 0;
			gridLine.X2 = width;
			gridLine.Y1 = gridLine.Y2 = height / 2;
			drawCanvas.Children.Add(gridLine);

			Line referenceLine;
			referenceLine = new Line() { Stroke = Brushes.Black, StrokeThickness = 3 };
			referenceLine.X1 = 10;
			referenceLine.X2 = 10 + 0.01 * xcoef;
			referenceLine.Y1 = referenceLine.Y2 = 10;
			drawCanvas.Children.Add(referenceLine);

			referenceLine = new Line() { Stroke = Brushes.Black, StrokeThickness = 3 };
			referenceLine.X1 = referenceLine.X2 = 10;
			referenceLine.Y1 = 10;
			referenceLine.Y2 = 10 + 0.01 * ycoef;
			drawCanvas.Children.Add(referenceLine);

			referenceLine = new Line() { Stroke = Brushes.Black, StrokeThickness = 3 };
			referenceLine.X1 = 20;
			referenceLine.X2 = 20 + 0.001 * xcoef;
			referenceLine.Y1 = referenceLine.Y2 = 20;
			drawCanvas.Children.Add(referenceLine);

			referenceLine = new Line() { Stroke = Brushes.Black, StrokeThickness = 3 };
			referenceLine.X1 = referenceLine.X2 = 20;
			referenceLine.Y1 = 20;
			referenceLine.Y2 = 20 + 0.001 * ycoef;
			drawCanvas.Children.Add(referenceLine);

			for (int i = 0; i < param.Lenss.Count; i++) {
				var Lens = param.Lenss[i];
				Line line = new Line() { Stroke = Brushes.Blue, StrokeThickness = 1 };
				line.X1 = Lens.distance * xcoef + xoffset;
				line.X2 = Lens.distance * xcoef + xoffset;
				line.Y1 = Lens.radius * ycoef + height / 2;
				line.Y2 = -Lens.radius * ycoef + height / 2;
				drawCanvas.Children.Add(line);

				double centerx = Lens.distance
					+ (Lens.isCurvatureRight ? 1 : -1) * Lens.thickness
					- (Lens.isCurvatureRight ^ Lens.curvature < 0 ? 1 : -1) * Math.Sqrt(Math.Pow(Lens.curvature, 2) - Math.Pow(Lens.radius, 2));

				const int maxCir = 360;
				bool drawEdge = false;
				for (int cir = -maxCir / 2; cir < maxCir / 2; cir++) {

					double x1, x2, y1, y2;
					double radOffset = (Lens.isCurvatureRight ? 0 : Math.PI);
					x1 = centerx + Lens.curvature * Math.Cos(cir * Math.PI / maxCir + radOffset);
					x2 = (centerx + Lens.curvature * Math.Cos((cir + 1) * Math.PI / maxCir + radOffset));
					y1 = Lens.curvature * Math.Sin(cir * Math.PI / maxCir + radOffset);
					y2 = Lens.curvature * Math.Sin((cir + 1) * Math.PI / maxCir + radOffset);

					if (Math.Abs(y1) > Lens.radius && Math.Abs(y2) < Lens.radius || (!drawEdge && Math.Abs(y1)<Lens.radius && Math.Abs(y2)<Lens.radius)) {
						Line lineu = new Line() { Stroke = Brushes.Blue, StrokeThickness = 1 };
						Line lineb = new Line() { Stroke = Brushes.Blue, StrokeThickness = 1 };
						lineu.X1 = Lens.distance * xcoef + xoffset;
						lineu.X2 = x2 * xcoef + xoffset;
						lineu.Y1 = Lens.radius * ycoef + height / 2;
						lineu.Y2 = Lens.radius * ycoef + height / 2;
						lineb.X1 = Lens.distance * xcoef + xoffset;
						lineb.X2 = x2 * xcoef + xoffset;
						lineb.Y1 = -Lens.radius * ycoef + height / 2;
						lineb.Y2 = -Lens.radius * ycoef + height / 2;
						drawCanvas.Children.Add(lineu);
						drawCanvas.Children.Add(lineb);
						drawEdge = true;
					}
					if (Math.Abs(y1) < Lens.radius && Math.Abs(y2) < Lens.radius) {
						Line linec = new Line() { Stroke = Brushes.Blue, StrokeThickness = 1 };
						linec.X1 = x1 * xcoef + xoffset;
						linec.X2 = x2 * xcoef + xoffset;
						linec.Y1 = y1 * ycoef + height / 2;
						linec.Y2 = y2 * ycoef + height / 2;
						drawCanvas.Children.Add(linec);
					}
				}
			}
		}
		private void paintResult(ParameterSet param, ResultSet res,Color col) {
			double width = drawCanvas.ActualWidth;
			double height = drawCanvas.ActualHeight;

			double xcoef = width / (maxDistance - minDistance);
			double ycoef = height / drawHeight;
			double xoffset = -minDistance * xcoef;


			// raytrace
			for (int i = 0; i < res.raytrace.Count; i++) {
				for (int j = 0; j + 1 < res.raytrace[i].Count; j++) {
					var pos = res.raytrace[i][j];
					var next = res.raytrace[i][j + 1];
					if (double.IsNaN(next.x) || double.IsNaN(next.y) || double.IsNaN(pos.x) || double.IsNaN(pos.y)) continue;
					//						Line line = new Line() { Stroke = Brushes.Black , StrokeThickness = 1 };
					Line line = new Line() {StrokeThickness = 1, Stroke = new SolidColorBrush(new Color() { 
						R = (byte)(col.R+(0x10 * i / res.raytrace.Count)), 
						G = (byte)(col.G+(0x10 - 0x10 * i / res.raytrace.Count)), 
						B = col.B, 
						A = 0xff 
					}) };
					line.X1 = pos.x * xcoef + xoffset;
					line.X2 = next.x * xcoef + xoffset;
					line.Y1 = -pos.y * ycoef + height / 2;
					line.Y2 = -next.y * ycoef + height / 2;
					drawCanvas.Children.Add(line);
				}
			}

			// wave front
			int wavefrontCount = 0;
			for (int i = 0; i < res.wavefront.Count; i++) {
				wavefrontCount = Math.Max(wavefrontCount, res.wavefront[i].Count);
			}
			for (int i = 0; i < wavefrontCount; i++) {
				for (int j = 0; j + 1 < res.wavefront.Count; j++) {
					if (res.wavefront[j].Count <= i || res.wavefront[j + 1].Count <= i) continue;
					var pos = res.wavefront[j][i];
					var next = res.wavefront[j + 1][i];
					if (double.IsNaN(next.x) || double.IsNaN(next.y) || double.IsNaN(pos.x) || double.IsNaN(pos.y)) continue;
					Line line = new Line() { Stroke = Brushes.Gray, StrokeThickness = 1 };
					line.X1 = pos.x * xcoef + xoffset;
					line.X2 = next.x * xcoef + xoffset;
					line.Y1 = -pos.y * ycoef + height / 2;
					line.Y2 = -next.y * ycoef + height / 2;
					drawCanvas.Children.Add(line);
				}
			}
		}

		private void Callback_Optmize(object sender, MouseButtonEventArgs e) {
/*
			DependencyObject dep = (DependencyObject)e.OriginalSource;

			// iteratively traverse the visual tree
			while ((dep != null) &&
					!(dep is DataGridCell) &&
					!(dep is DataGridColumnHeader))
			{
				dep = VisualTreeHelper.GetParent(dep);
			}

			if (dep == null)
				return;

			if (dep is DataGridColumnHeader){
				DataGridColumnHeader columnHeader = dep as DataGridColumnHeader;
				// do something
			}

			if (dep is DataGridCell) {
				DataGridCell cell = dep as DataGridCell;
				while ((dep != null) && !(dep is DataGridRow)) {
					dep = VisualTreeHelper.GetParent(dep);
				}
				DataGridRow row = dep as DataGridRow;

				DataGrid dataGrid = (DataGrid)ItemsControl.ItemsControlFromItemContainer(row);
				int rowIndex = dataGrid.ItemContainerGenerator.IndexFromContainer(row);
				int columnIndex = cell.Column.DisplayIndex;

				if (rowIndex < LensData.Count) {

					double org = double.Parse(LensData[rowIndex].distance);

					List<double> coefs = new List<double>() { 1.0, 1.1, 0.9, 1.01, 0.99, 1.001, 0.999 };
					int bestCoef = 0;
					double bestAbberation = double.MaxValue;

					for (int i = 0; i < coefs.Count; i++) {
						LensData[rowIndex].distance = (org * coefs[i]).ToString();
						createParameterSet();
						res = (new LensAnalyze()).solve(param);
//						if (res.lastAbberationWidth < bestAbberation) {
//							bestAbberation = res.lastAbberationWidth;
//							bestCoef = i;
//						}
						if (res.lastAbberationHeight < bestAbberation) {
							bestAbberation = res.lastAbberationHeight;
							bestCoef = i;
						}
					}
					LensData[rowIndex].distance = (org * coefs[bestCoef]).ToString();
					if (bestCoef != 0) update();
				}
			}
*/
		}
	}
}
