import { useState, ChangeEvent } from "react";
import "../App.css";

const SubmitJobs: React.FC = () => {
  const [email, setEmail] = useState("");
  const [jobName, setJobName] = useState("");
  const [pdbFile, setPdbFile] = useState<File | null>(null);
  const [submissionSuccess, setSubmissionSuccess] = useState(false);

  const handlePdbFileChange = (e: ChangeEvent<HTMLInputElement>) => {
    if (e.target.files && e.target.files.length > 0) {
      setPdbFile(e.target.files[0]);
    }
  };

  const fetchDefaultFile = async () => {
    console.log("Fetching default file...");
    const response = await fetch("https://files.rcsb.org/download/2KL8.pdb");

    const blob = await response.blob();
    const file = new File([blob], "2KL8.pdb", {
      type: "application/octet-stream",
    });
    setPdbFile(file);
  };

  const handleJobSubmit = async (e: { preventDefault: () => void }) => {
    e.preventDefault();
    if (!email || !jobName || !pdbFile) {
      alert("Please provide both a job name and a PDB file.");
      return;
    }

    const toBase64 = (file: Blob) =>
      new Promise((resolve, reject) => {
        const reader = new FileReader();
        reader.readAsDataURL(file);
        reader.onload = () => {
          if (reader.result) {
            resolve((reader.result as string).split(",")[1]);
          } else {
            reject(new Error("File reading failed"));
          }
        };
        reader.onerror = (error) => reject(error);
      });

    const fileBase64 = await toBase64(pdbFile);

    const payload = {
      email,
      jobName,
      file: fileBase64,
    };

    try {
      const response = await fetch(
        "https://w0kfu9l8nk.execute-api.us-east-1.amazonaws.com/default/submitJob",
        {
          method: "POST",
          headers: { "Content-Type": "application/json" },
          body: JSON.stringify(payload),
        }
      );

      const data = await response.json();
      if (response.ok) {
        setSubmissionSuccess(true); // Show success animation
        console.log("Job submitted successfully!");

        // Reset form fields after animation delay
        setTimeout(() => {
          setJobName("");
          setPdbFile(null);
          setSubmissionSuccess(false);
        }, 3000); // 3 seconds for the animation to display
      } else {
        alert("Failed to submit job: " + data.error);
      }
    } catch (error) {
      if (error instanceof Error) {
        alert("An error occurred: " + error.message);
      } else {
        alert("An unknown error occurred");
      }
    }
  };

  return (
    <div className="right-side">
      {/* New Section for Job Name and PDB File Upload */}
      <h2>Prototype RFdiffusion + Gromacs automation with Zarrtraj</h2>
      <div className="instructions">
        <p>
          This tool will run RFdiffusion on your protein to attempt to generate
          a similar protein with different residues and folding. Then, it will
          solvate the resulting protein in a water box and run a short MD
          simulation to observe its behavior. The results will be provided to
          you via S3 URLs on the next page once the job is complete. Topologies
          and trajectories can be streamed directly from the S3 URL using{" "}
          <a href="https://github.com/becksteinlab/zarrtraj">Zarrtraj</a>. See
          the final page for more.
        </p>
      </div>
      <input
        type="text"
        placeholder="Job Name"
        value={jobName}
        onChange={(e) => setJobName(e.target.value)}
      />
      <div className="file-upload">
        <div>
          <input
            id="pdb-submit"
            type="file"
            accept=".pdb"
            onChange={handlePdbFileChange}
            style={{ display: "none" }}
          />
          <label htmlFor="pdb-submit" className="default-file-button">
            Select PDB File
          </label>
        </div>
        <button className="default-file-button" onClick={fetchDefaultFile}>
          Use sample file
        </button>
        {pdbFile && <span className="fname">{pdbFile.name}</span>}
      </div>
      <input
        type="text"
        placeholder="Email address"
        value={email}
        onChange={(e) => setEmail(e.target.value)}
      />
      <button className="submit-button" onClick={handleJobSubmit}>
        Submit Job
      </button>

      {/* Green check animation */}
      {submissionSuccess && (
        <div className="success-animation">
          <div className="checkmark"></div>
        </div>
      )}
    </div>
  );
};

export default SubmitJobs;
