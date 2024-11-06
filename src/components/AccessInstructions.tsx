import { useState, ChangeEvent } from "react";
import "../App.css";
import { CopyBlock, dracula } from "react-code-blocks";

const AccessInstructions: React.FC = () => {
  const installationText = `conda create -c conda-forge -n zarrtraj zarrtraj MDAnalysis`;
  const codeText = "";

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
        and
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
