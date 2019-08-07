import org.jlab.groot.data.TDirectory
import org.jlab.groot.group.DataGroup;

public class Filemanager{
	public void main(){
		this.dir_in = new TDirectory ();
		this.dir_out= new TDirectory ();
	}
	public void readhipo(String filename){
		this.dir_in.readfile(filename)
	}
	public void getObject(String objectname){
		this.dir_in.getObject(objectname)
	}
	public void mkdir(args){
		this.dir
	}
	public void writehipo(String filename){
		this.dir_out.writefile(filename)
	}
	public void addDataSet(args){
		for arg in args{
			this.dir.addDataSet(arg)
		}
	}
}

