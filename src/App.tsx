import { useState, ChangeEvent } from 'react'
import './App.css'

function App() {
  const [email, setEmail] = useState('');
  const [jobName, setJobName] = useState('');
  const [pdbFile, setPdbFile] = useState<File | null>(null);

  const handlePdbFileChange = (e: ChangeEvent<HTMLInputElement>) => {
    if (e.target.files && e.target.files.length > 0) {
      setPdbFile(e.target.files[0]);
    }
  };

  const handleJobSubmit = () => {
    if (!jobName || !pdbFile) {
      alert('Please provide both a job name and a PDB file.');
      return;
    }

    // Process job submission logic here (e.g., upload to server)
    console.log('Job Name:', jobName);
    console.log('PDB File:', pdbFile);

    // Reset form fields
    setJobName('');
    setPdbFile(null);
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
          <input type="file" accept=".pdb" onChange={handlePdbFileChange} />
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
  )
}

export default App
