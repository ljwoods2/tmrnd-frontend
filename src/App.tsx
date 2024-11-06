import { useState, ChangeEvent, useEffect } from "react";
import "./App.css";

function App() {
  const [email, setEmail] = useState("");
  const [jobName, setJobName] = useState("");
  const [pdbFile, setPdbFile] = useState<File | null>(null);

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

  const handleJobSubmit = async (e) => {
    e.preventDefault();
    if (!email || !jobName || !pdbFile) {
      alert("Please provide both a job name and a PDB file.");
      return;
    }

    const toBase64 = (file) =>
      new Promise((resolve, reject) => {
        const reader = new FileReader();
        reader.readAsDataURL(file);
        reader.onload = () => resolve(reader.result.split(",")[1]); // Remove the base64 prefix
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
        alert("Job submitted successfully!");
        console.log("Job Name:", jobName);
        console.log("PDB File:", pdbFile);

        // Reset form fields
        setJobName("");
        setPdbFile(null);
      } else {
        alert("Failed to submit job: " + data.error);
      }
    } catch (error) {
      alert("An error occurred: " + error.message);
    }
  };

  return (
    <div className="background-gradient">
      <div className="container">
        <div className="right-side">
          {/* New Section for Job Name and PDB File Upload */}
          <h2>Submit Job</h2>
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
        </div>
      </div>
    </div>
  );
}

export default App;
