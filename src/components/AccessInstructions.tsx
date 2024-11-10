import "../App.css";
import { CopyBlock, dracula } from "react-code-blocks";

const AccessInstructions: React.FC = () => {
  const installationText = `conda create -c conda-forge -n zarrtraj MDAnalysis`;
  const codeText = `import zarrtraj\nimport MDAnalysis as mda\nimport fsspec\nimport MDAnalysis.analysis.rms
with fsspec.open("s3://tmrnd-prototype/job_name_diffused.pdb",
  \t\t\t"r", anon=True) as top:
  \tstorage_options = {"anon": True}
  \tu = mda.Universe(top,
  \t\t"s3://tmrnd-prototype/job_name_diffused.zarrmd",
  \t\tstorage_options = {"anon": True},
  \t\ttopology_format="PDB")
  \tR = MDAnalysis.analysis.rms.RMSD(u, u,
           select="backbone")
  \tR.run()
  \tresults = R.rmsd`;

  return (
    <div className="usage">
      <h2>Installation Instructions</h2>
      <div className="instructions">
        <div>
          <p>
            This prototype uses{" "}
            <a href="https://github.com/becksteinlab/zarrtraj">Zarrtraj</a> to
            serve trajectory data. To access the data, you will need to install
            the Zarrtraj package and use the provided Python code to load the
            trajectory data.
          </p>
        </div>
        <div className="code-block">
          <CopyBlock
            text={installationText}
            language="bash"
            showLineNumbers={true}
            theme={dracula}
          />
        </div>
        Once installed, you can interact with the trajectory data like you
        typically would with MDAnalysis. In this example, we load a trajectory
        and perform RMSD on the backbone of the protein.
        <div className="code-block">
          <CopyBlock
            text={codeText}
            language="python"
            showLineNumbers={true}
            theme={dracula}
          />
        </div>
      </div>
    </div>
  );
};

export default AccessInstructions;
