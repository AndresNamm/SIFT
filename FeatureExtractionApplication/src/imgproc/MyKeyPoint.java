package imgproc;

import org.opencv.features2d.KeyPoint;

import java.util.ArrayList;
import java.util.List;

import static java.lang.Math.pow;

/**
 * Created by UDU on 31.07.2016.
 */
public class MyKeyPoint {

	public int getY() {
		return y;
	}

	public int getX() {
		return x;
	}

	public int getOctave() {
		return octave;
	}

	public int getScale() {
		return scale;
	}

	public int getLayer() {
		return layer;
	}

	public double getSigma() {
		return sigma;
	}

	public float getSize() {
		return size;
	}

	public float getResponse() {
		return response;
	}

	private List<Float> angle_deg = new ArrayList<>();
	private List<Double> angle_mag = new ArrayList<>();

	private int y = 0;// height
	private int x = 0; // size
	private int octave = 0;
	private int scale = 0;
	private int layer = 0;
	private double sigma;

	private float size = 2;
	private float angle = 0;
	private float response = 9999;

	public MyKeyPoint(int x, int y, int octave, int scale, double response) {

		this.x = (int) (x * pow(2, octave));
		this.y = (int) (y * pow(2, octave));
		this.octave = octave;
		this.response = (float) response;
		this.layer = scale;
		// this.size = sigma * Math.pow(2.0, (layer)) / nOctaveLayers) * (1 <<
		// octv) * 2;

	}

	public boolean setAngle(float d) {
		this.angle = d;
		return true;
	}

	public boolean setSigma(float d) {
		this.sigma = d;
		return true;
	}

	public boolean addSize(float d) {
		this.size = d;
		System.out.println("Current Sigma" + sigma);
		// this.size = (float) Math.abs((sigma * Math.pow(2.0, (layer)) / 2.0) *
		// (1 << octave) * 2 ) ;
		return true;
	}

	public KeyPoint returnSimple() {
		return new KeyPoint(x, y, size);
	}

	public int getListsize() {
		return angle_deg.size();
	}

	public ArrayList<KeyPoint> getOpenCvKp() {
		ArrayList<KeyPoint> temp = new <KeyPoint>ArrayList();

		for (int i = 0; i < angle_deg.size(); i++) {
			temp.add(new KeyPoint(this.x, this.y, this.size, this.angle_deg.get(i), this.response, this.octave));
			// System.out.println(this.response);
		}
		return temp;
	}

	public void addAngle(float angle) {
		angle_deg.add(angle);
	}

	public void addAngleMag(Double mag) {
		angle_mag.add(mag);
	}

}
